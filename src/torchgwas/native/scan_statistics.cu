// Optional shared OLS elementwise/reduction kernels; GEMM remains in Torch.
#include <cuda_runtime.h>
#include <cstdint>
#include <cmath>
#include <string>
namespace {
thread_local std::string last_error;
constexpr unsigned threads=256;
template<typename T> __device__ float dosage(T value,float scale,int missing,float sentinel){
  float x=static_cast<float>(value);
  return missing && x==sentinel ? nanf("") : x/scale;
}
// Which two-bit code table a packed row uses. PGEN stores the dosage
// directly with 3 meaning missing; PLINK1 orders its codes differently
// and reserves 1 for missing.
enum PackedCodes { PACKED_NONE = 0, PACKED_PGEN = 1, PACKED_PLINK1 = 2 };
// THE TWO BRANCHES COUNT DIFFERENT ALLELES, AND EACH IS CORRECT FOR ITS FORMAT.
// PLINK1 decodes to A2 dosage (.bim column 6); PGEN decodes to ALT. Verified
// against the official PLINK example (cog-genomics .bed spec): bim "G A", PED
// GG/AA/missing/AA/AA/AA decodes to [0, 2, nan, 2, 2, 2] and the reader reports
// effect allele "A" -- self-consistent, counting and declaring the same allele.
//
// For a .bim written from a .pvar, A2 == REF and ALT == A1, so the SAME cohort
// read as .bed and as .pgen yields dosages summing to exactly 2.0 per variant
// and betas of opposite sign. That is a CONTRACT DIFFERENCE, not a decode bug:
// each reader names its own effect allele. Making them agree means flipping the
// dosage AND the declared effect allele together -- flipping only the dosage
// (attempted 2026-09-15) leaves the reader counting A1 while announcing A2,
// which is worse than either consistent choice.
__device__ __forceinline__ float decode_two_bit(unsigned code, int kind){
  if(kind==PACKED_PLINK1) return code==1u?nanf(""):float((code+1u)>>1);
  return code==3u?nanf(""):float(code);
}
template<typename T, int packed=PACKED_NONE> __global__ void prepare_kernel(const T* input,int64_t samples,float scale,
    int missing,float sentinel,float* centered,float* ss,float* minimum,float* maximum,
    int32_t* present_out,int64_t row_stride=0){
  const int64_t variant=blockIdx.x;
  const T* row=input+variant*(packed!=PACKED_NONE?row_stride:samples);
  __shared__ double sums[threads];
  __shared__ float lows[threads],highs[threads];
  __shared__ unsigned counts[threads];
  double sum=0;float lo=INFINITY,hi=-INFINITY;unsigned observed=0;
  for(int64_t sample=threadIdx.x;sample<samples;sample+=blockDim.x){
    float value;
    if constexpr(packed){unsigned code=(unsigned(row[sample/4])>>(2*(sample%4)))&3u;value=decode_two_bit(code,packed);}
    else value=dosage(row[sample],scale,missing,sentinel);
    // A missing call takes no part in the mean or the range. Counting the
    // observed calls here costs nothing: this pass already reads every value,
    // and the count is what the second pass needs to centre correctly.
    if(!isnan(value)){sum+=double(value);lo=fminf(lo,value);hi=fmaxf(hi,value);++observed;}
  }
  sums[threadIdx.x]=sum;lows[threadIdx.x]=lo;highs[threadIdx.x]=hi;counts[threadIdx.x]=observed;
  __syncthreads();
  for(unsigned stride=threads/2;stride;stride>>=1){
    if(threadIdx.x<stride){sums[threadIdx.x]+=sums[threadIdx.x+stride];
      lows[threadIdx.x]=fminf(lows[threadIdx.x],lows[threadIdx.x+stride]);
      highs[threadIdx.x]=fmaxf(highs[threadIdx.x],highs[threadIdx.x+stride]);
      counts[threadIdx.x]+=counts[threadIdx.x+stride];}
    __syncthreads();
  }
  const unsigned present=counts[0];
  const float mean=present?static_cast<float>(sums[0]/double(present)):0.0f;
  // With nothing observed the range stays empty, so maximum>minimum fails
  // downstream and the variant is reported invariant rather than silently
  // producing a statistic from no data.
  // The count the mean was taken over is what the variant's df is made of.
  if(threadIdx.x==0){minimum[variant]=present?lows[0]:nanf("");maximum[variant]=present?highs[0]:nanf("");
    present_out[variant]=static_cast<int32_t>(present);}
  double squares=0;
  for(int64_t sample=threadIdx.x;sample<samples;sample+=blockDim.x){
    float value;
    if constexpr(packed){unsigned code=(unsigned(row[sample/4])>>(2*(sample%4)))&3u;value=decode_two_bit(code,packed);}
    else value=dosage(row[sample],scale,missing,sentinel);
    // Centre, and mask: a missing call becomes exactly zero, so it drops out
    // of the GEMM and the sum of squares without a separate mask tensor.
    value=isnan(value)?0.0f:value-mean;
    centered[variant*samples+sample]=value;squares+=double(value)*double(value);
  }
  // All threads have consumed sums[0] before it is reused for the SS reduction.
  __syncthreads();sums[threadIdx.x]=squares;__syncthreads();
  for(unsigned stride=threads/2;stride;stride>>=1){
    if(threadIdx.x<stride)sums[threadIdx.x]+=sums[threadIdx.x+stride];__syncthreads();
  }
  if(threadIdx.x==0)ss[variant]=static_cast<float>(sums[0]);
}
__global__ void finish_kernel(const float* products,const float* centered_ss,const float* minimum,
    const float* maximum,const float* phenotype_ss,const int32_t* present,float df_offset,
    int64_t traits,int64_t covariates,
    int validate_range,float* beta,float* statistic,uint8_t* status){
  const int64_t variant=blockIdx.x;
  const float* row=products+variant*(traits+covariates);
  __shared__ float sums[threads];
  float sum=0;
  for(int64_t c=threadIdx.x;c<covariates;c+=blockDim.x){float x=row[traits+c];sum+=x*x;}
  sums[threadIdx.x]=sum;__syncthreads();
  for(unsigned stride=threads/2;stride;stride>>=1){
    if(threadIdx.x<stride)sums[threadIdx.x]+=sums[threadIdx.x+stride];__syncthreads();
  }
  const float residual=centered_ss[variant]-sums[0];
  // Each variant spends its own observed samples, so it keeps its own
  // residual degrees of freedom. One with none left is not a result.
  const float df=static_cast<float>(present[variant])+df_offset;
  const bool valid=residual>1e-12f && maximum[variant]>minimum[variant] && df>0.0f;
  if(threadIdx.x==0){
    uint8_t code=!isfinite(residual)?1:(valid?0:2);
    if(validate_range && (minimum[variant]<0 || maximum[variant]>2))code=3;
    status[variant]=code;
  }
  const float safe=isnan(residual)?residual:fmaxf(residual,1e-12f);
  for(int64_t trait=threadIdx.x;trait<traits;trait+=blockDim.x){
    float gy=row[trait],b=gy/safe;
    float yss=phenotype_ss[trait]-gy*gy/safe;
    yss=isnan(yss)?yss:fmaxf(yss,1e-12f);
    float se=sqrtf(yss/df/safe);
    beta[variant*traits+trait]=valid?b:0;
    statistic[variant*traits+trait]=valid?b/se:0;
  }
}
int checked_launch(){auto error=cudaGetLastError();if(error!=cudaSuccess){last_error=cudaGetErrorString(error);return -1;}return 0;}
}
extern "C" {
int tg_scan_abi_version(){return 4;}
const char* tg_scan_error(){return last_error.c_str();}
int tg_scan_prepare(const void* input,int dtype,int64_t variants,int64_t samples,float scale,
    int missing,float sentinel,float* centered,float* ss,float* minimum,float* maximum,
    int32_t* present,void* stream){
  if(variants<=0 || samples<=0 || !std::isfinite(scale) || scale<=0){last_error="invalid preparation dimensions/scale";return -1;}
  cudaStream_t s=static_cast<cudaStream_t>(stream);
  if(dtype==0)prepare_kernel<<<variants,threads,0,s>>>(static_cast<const int8_t*>(input),samples,scale,missing,sentinel,centered,ss,minimum,maximum,present);
  else if(dtype==1)prepare_kernel<<<variants,threads,0,s>>>(static_cast<const uint8_t*>(input),samples,scale,missing,sentinel,centered,ss,minimum,maximum,present);
  else if(dtype==2)prepare_kernel<<<variants,threads,0,s>>>(static_cast<const float*>(input),samples,scale,missing,sentinel,centered,ss,minimum,maximum,present);
  else{last_error="unsupported preparation dtype";return -1;}
  return checked_launch();
}
int tg_scan_prepare_pgen2(const uint8_t* input,int64_t variants,int64_t samples,int64_t row_stride,
    float* centered,float* ss,float* minimum,float* maximum,int32_t* present,void* stream){
  if(variants<=0 || samples<=0 || row_stride<(samples+3)/4){last_error="invalid packed PGEN dimensions";return -1;}
  prepare_kernel<uint8_t,PACKED_PGEN><<<variants,threads,0,static_cast<cudaStream_t>(stream)>>>(input,samples,1,0,0,centered,ss,minimum,maximum,present,row_stride);
  return checked_launch();
}
// Same two-bit layout, different code table: PLINK1 reserves 01 for missing
// and orders the rest so the dosage is (code + 1) >> 1.
int tg_scan_prepare_bed2(const uint8_t* input,int64_t variants,int64_t samples,int64_t row_stride,
    float* centered,float* ss,float* minimum,float* maximum,int32_t* present,void* stream){
  if(variants<=0 || samples<=0 || row_stride<(samples+3)/4){last_error="invalid packed BED dimensions";return -1;}
  prepare_kernel<uint8_t,PACKED_PLINK1><<<variants,threads,0,static_cast<cudaStream_t>(stream)>>>(input,samples,1,0,0,centered,ss,minimum,maximum,present,row_stride);
  return checked_launch();
}
int tg_scan_finish(const float* products,const float* ss,const float* minimum,const float* maximum,
    const float* phenotype_ss,const int32_t* present,float df_offset,
    int64_t variants,int64_t traits,int64_t covariates,
    int validate_range,float* beta,float* statistic,uint8_t* status,void* stream){
  if(variants<=0 || traits<=0 || covariates<0 || !std::isfinite(df_offset)){last_error="invalid statistics dimensions/df";return -1;}
  finish_kernel<<<variants,threads,0,static_cast<cudaStream_t>(stream)>>>(products,ss,minimum,maximum,phenotype_ss,present,df_offset,traits,covariates,validate_range,beta,statistic,status);
  return checked_launch();
}
}
