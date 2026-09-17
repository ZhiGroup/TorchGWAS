// Direct device-resident Layout-2 decoder, derived from the validated converter.
// CPU batch preparation is independent of GPU decode and contains no store writer.
#include <cuda_runtime.h>
#include <nvcomp/deflate.h>
#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fcntl.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
namespace {

constexpr uint32_t kAdlerMod = 65521;

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }

void check_cuda(cudaError_t status, const char* what) {
  if (status != cudaSuccess) fail(std::string(what) + ": " + cudaGetErrorString(status));
}

void check_nvcomp(nvcompStatus_t status, const char* what) {
  if (status != nvcompSuccess) fail(std::string(what) + ": " + nvcompGetStatusString(status));
}

uint16_t le16(const uint8_t* p) { return uint16_t(p[0]) | (uint16_t(p[1]) << 8); }
uint32_t le32(const uint8_t* p) {
  return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}
uint32_t be32(const uint8_t* p) {
  return (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | uint32_t(p[3]);
}
size_t align_up(size_t x, size_t a) { return (x + a - 1) / a * a; }

void pread_full(int fd, void* buffer, size_t bytes, uint64_t offset) {
  auto* p = static_cast<uint8_t*>(buffer);
  size_t done = 0;
  while (done < bytes) {
    ssize_t got = ::pread(fd, p + done, bytes - done, static_cast<off_t>(offset + done));
    if (got < 0) fail("pread failed at byte " + std::to_string(offset + done) + ": " + std::strerror(errno));
    if (got == 0) fail("short BGEN read at byte " + std::to_string(offset + done));
    done += static_cast<size_t>(got);
  }
}

struct Options {
  size_t n_bgen_samples = 0, n_samples = 0;
  int gpu = 0;
  // Verifying every record's stored Adler-32 costs about an eighth of decode.
  // On by default: it is the only end-to-end integrity check on this path.
  bool verify_adler = true;
};

struct GpuRecord {
  size_t comp_offset = 0, comp_size = 0, raw_offset = 0, raw_size = 0;
  uint32_t expected_adler = 0;
};

struct Batch {
  std::vector<uint8_t> compressed;
  std::vector<GpuRecord> records;
};

__device__ uint32_t read_le32_device(const uint8_t* p) {
  return uint32_t(p[0]) | (uint32_t(p[1]) << 8) | (uint32_t(p[2]) << 16) | (uint32_t(p[3]) << 24);
}
__device__ uint16_t read_le16_device(const uint8_t* p) { return uint16_t(p[0]) | (uint16_t(p[1]) << 8); }

__device__ uint64_t read_bits(const uint8_t* packed, uint64_t bit_offset, int bits) {
  const uint64_t byte_offset = bit_offset >> 3;
  const int shift = int(bit_offset & 7);
  const int bytes = (shift + bits + 7) >> 3;
  uint64_t word = 0;
  #pragma unroll
  for (int j = 0; j < 5; ++j) if (j < bytes) word |= uint64_t(packed[byte_offset + j]) << (8 * j);
  const uint64_t mask = bits == 32 ? 0xffffffffULL : ((1ULL << bits) - 1ULL);
  return (word >> shift) & mask;
}

__global__ void decode_layout2_kernel(
    const uint8_t* const* raw_ptrs, const size_t* actual_sizes, const size_t* expected_sizes,
    const int* nv_status, const int32_t* sample_indices, size_t n_selected,
    size_t n_bgen_samples, float* output, int* parse_status, size_t n_records) {
  const size_t v = blockIdx.x;
  if (v >= n_records) return;
  const uint8_t* raw = raw_ptrs[v];
  if (threadIdx.x == 0) {
    int status = 0;
    if (nv_status[v] != int(nvcompSuccess)) status = 1;
    else if (actual_sizes[v] != expected_sizes[v] || actual_sizes[v] < 10 + n_bgen_samples) status = 2;
    else {
      const uint32_t n = read_le32_device(raw);
      const uint16_t k = read_le16_device(raw + 4);
      const uint8_t min_ploidy = raw[6], max_ploidy = raw[7];
      const uint8_t phased = raw[8 + n_bgen_samples];
      const uint8_t bits = raw[9 + n_bgen_samples];
      if (n != n_bgen_samples || k != 2) status = 3;
      else if (min_ploidy != 2 || max_ploidy != 2) status = 4;
      else if (phased != 0) status = 6;
      else if (bits < 1 || bits > 32) status = 7;
      else {
        const uint64_t packed_bytes = (uint64_t(n_bgen_samples) * 2 * bits + 7) / 8;
        if (10 + n_bgen_samples + packed_bytes != actual_sizes[v]) status = 8;
      }
    }
    parse_status[v] = status;
  }
  __syncthreads();
  if (parse_status[v] != 0) return;
  const int bits = raw[9 + n_bgen_samples];
  const uint64_t denom = bits == 32 ? 0xffffffffULL : ((1ULL << bits) - 1ULL);
  const uint8_t* packed = raw + 10 + n_bgen_samples;
  for (size_t j = threadIdx.x; j < n_selected; j += blockDim.x) {
    const uint32_t s = static_cast<uint32_t>(sample_indices[j]);
    const uint8_t pm = raw[8 + s];
    if ((pm & 63) != 2) { atomicCAS(parse_status + v, 0, 4); continue; }
    if (pm & 128) { output[v*n_selected+j] = __int_as_float(0x7fffffff); continue; }
    const uint64_t p0 = read_bits(packed, uint64_t(s) * 2 * bits, bits);
    const uint64_t p1 = read_bits(packed, uint64_t(s) * 2 * bits + bits, bits);
    if (p0 + p1 > denom) { atomicCAS(parse_status + v, 0, 9); continue; }
    const uint64_t dosage_num = 2 * denom - 2 * p0 - p1;
    output[v * n_selected + j] = static_cast<float>(double(dosage_num) / double(denom));
  }
}

__global__ void adler32_kernel(
    const uint8_t* const* raw_ptrs, const size_t* sizes, const uint32_t* expected,
    const int* nv_status, int* parse_status, size_t n_records) {
  const size_t v = blockIdx.x;
  if (v >= n_records || nv_status[v] != int(nvcompSuccess) || parse_status[v] != 0) return;
  const uint8_t* data = raw_ptrs[v];
  const uint64_t n = sizes[v];
  // Adler-32 is A = 1 + sum(b_i), B = n + sum((n-i)*b_i). The weight depends
  // only on the global index, so threads may take any partition of the bytes
  // and still compute their exact contribution. Striding by blockDim keeps
  // consecutive threads on consecutive bytes: one contiguous 256-byte span
  // per iteration instead of 256 separate positions megabytes apart, which
  // is what held this kernel to roughly 2% of achievable bandwidth.
  //
  // No intermediate reduction is needed: for a record of n bytes the
  // weighted sum is bounded by 255*n*(n+1)/2, which stays inside uint64 for
  // any BGEN record (n is capped by the decode workspace, far below the
  // ~1.2e8 bytes where this would overflow).
  uint64_t sum = 0, weighted = 0;
  for (uint64_t i = threadIdx.x; i < n; i += blockDim.x) {
    const uint64_t b = data[i];
    sum += b;
    weighted += (n - i) * b;
  }
  __shared__ uint64_t sums[256];
  __shared__ uint64_t weights[256];
  sums[threadIdx.x] = sum % kAdlerMod;
  weights[threadIdx.x] = weighted % kAdlerMod;
  __syncthreads();
  for (unsigned stride = blockDim.x / 2; stride; stride >>= 1) {
    if (threadIdx.x < stride) {
      sums[threadIdx.x] = (sums[threadIdx.x] + sums[threadIdx.x + stride]) % kAdlerMod;
      weights[threadIdx.x] = (weights[threadIdx.x] + weights[threadIdx.x + stride]) % kAdlerMod;
    }
    __syncthreads();
  }
  if (threadIdx.x == 0) {
    const uint32_t a = static_cast<uint32_t>((1 + sums[0]) % kAdlerMod);
    const uint32_t b = static_cast<uint32_t>((n + weights[0]) % kAdlerMod);
    const uint32_t observed = (b << 16) | a;
    if (observed != expected[v]) atomicCAS(parse_status + v, 0, 10);
  }
}

template <typename T>
void ensure_device(T*& ptr, size_t& capacity_items, size_t required_items) {
  if (capacity_items >= required_items) return;
  const size_t next_capacity=std::max(required_items, capacity_items + capacity_items / 2 + 1);
  T* old=ptr;ptr=nullptr;capacity_items=0;
  if (old) check_cuda(cudaFree(old), "cudaFree grow");
  check_cuda(cudaMalloc(reinterpret_cast<void**>(&ptr), next_capacity * sizeof(T)), "cudaMalloc grow");
  capacity_items=next_capacity;
}

template <typename T>
void ensure_pinned(T*& ptr, size_t& capacity_items, size_t required_items) {
  if (capacity_items >= required_items) return;
  const size_t next_capacity=std::max(required_items, capacity_items + capacity_items / 2 + 1);
  T* old=ptr;ptr=nullptr;capacity_items=0;
  if (old) check_cuda(cudaFreeHost(old), "cudaFreeHost grow");
  check_cuda(cudaMallocHost(reinterpret_cast<void**>(&ptr), next_capacity * sizeof(T)), "cudaMallocHost grow");
  capacity_items=next_capacity;
}

struct GpuBuffers {
  uint8_t *d_comp = nullptr, *d_raw = nullptr, *d_output = nullptr, *d_temp = nullptr;
  const void **d_comp_ptrs = nullptr;
  void **d_raw_ptrs = nullptr;
  size_t *d_comp_sizes = nullptr, *d_raw_caps = nullptr, *d_actual_sizes = nullptr;
  uint32_t* d_adler = nullptr;
  int *d_nv_status = nullptr, *d_parse_status = nullptr;
  int32_t* d_sample_indices = nullptr;
  uint8_t *h_comp = nullptr, *h_output = nullptr;
  int* h_parse_status = nullptr;
  size_t comp_cap=0, raw_cap=0, output_cap=0, temp_cap=0, record_cap=0, sample_cap=0;
  size_t h_comp_cap=0, h_output_cap=0, h_status_cap=0;
  cudaStream_t stream{};
  ~GpuBuffers() {
    cudaFree(d_comp); cudaFree(d_raw); cudaFree(d_output); cudaFree(d_temp);
    cudaFree(d_comp_ptrs); cudaFree(d_raw_ptrs); cudaFree(d_comp_sizes); cudaFree(d_raw_caps);
    cudaFree(d_actual_sizes); cudaFree(d_adler); cudaFree(d_nv_status); cudaFree(d_parse_status);
    cudaFree(d_sample_indices); cudaFreeHost(h_comp); cudaFreeHost(h_output); cudaFreeHost(h_parse_status);
    if (stream) cudaStreamDestroy(stream);
  }
};

void ensure_record_buffers(GpuBuffers& g, size_t n) {
  if (g.record_cap >= n) return;
  const size_t capacity=std::max(n, g.record_cap + g.record_cap / 2 + 1);
  g.record_cap=0;
  auto release=[](auto*& pointer){auto* old=pointer;pointer=nullptr;if(old)check_cuda(cudaFree(old),"free record buffer");};
  release(g.d_comp_ptrs);release(g.d_raw_ptrs);release(g.d_comp_sizes);release(g.d_raw_caps);
  release(g.d_actual_sizes);release(g.d_adler);release(g.d_nv_status);release(g.d_parse_status);
  check_cuda(cudaMalloc(&g.d_comp_ptrs, capacity * sizeof(void*)), "alloc comp ptrs");
  check_cuda(cudaMalloc(&g.d_raw_ptrs, capacity * sizeof(void*)), "alloc raw ptrs");
  check_cuda(cudaMalloc(&g.d_comp_sizes, capacity * sizeof(size_t)), "alloc comp sizes");
  check_cuda(cudaMalloc(&g.d_raw_caps, capacity * sizeof(size_t)), "alloc raw sizes");
  check_cuda(cudaMalloc(&g.d_actual_sizes, capacity * sizeof(size_t)), "alloc actual sizes");
  check_cuda(cudaMalloc(&g.d_adler, capacity * sizeof(uint32_t)), "alloc adler");
  check_cuda(cudaMalloc(&g.d_nv_status, capacity * sizeof(int)), "alloc nv statuses");
  check_cuda(cudaMalloc(&g.d_parse_status, capacity * sizeof(int)), "alloc parse statuses");
  g.record_cap=capacity;
}


// C ABI keeps the extension independent of the installed PyTorch C++ ABI.
// Caller owns output. This decoder only copies compact statuses back to host.
thread_local std::string last_error;
struct Decoder {
  Options options;
  GpuBuffers gpu;
  bool profile = false;
  cudaEvent_t timings[6]{};
  cudaEvent_t completion{};
  double totals[16]{};
  Decoder(int device, size_t n, const int32_t* selection, size_t selected) {
    options.gpu = device; options.n_bgen_samples = n; options.n_samples = selected;
    check_cuda(cudaSetDevice(device), "select decoder GPU");
    for (size_t i=0;i<selected;++i) if (selection[i]<0 || size_t(selection[i])>=n) fail("sample index out of bounds");
    ensure_device(gpu.d_sample_indices, gpu.sample_cap, selected);
    check_cuda(cudaMemcpy(gpu.d_sample_indices, selection, selected*sizeof(int32_t), cudaMemcpyHostToDevice), "sample selection");
    check_cuda(cudaStreamCreateWithFlags(&gpu.stream, cudaStreamNonBlocking), "decoder stream");
    const char* setting=std::getenv("TORCHGWAS_BGEN_PROFILE");
    profile=setting && std::strcmp(setting,"0")!=0;
    try {
    const char* blocking=std::getenv("TORCHGWAS_BGEN_BLOCKING_SYNC");
    if(blocking && std::strcmp(blocking,"1")==0)
      check_cuda(cudaEventCreateWithFlags(&completion,cudaEventBlockingSync|cudaEventDisableTiming),"blocking completion event");
    totals[0]=profile?1:0;
    if(profile) for(auto& event:timings) check_cuda(cudaEventCreate(&event),"profile event");
    } catch (...) {
      for(auto event:timings)if(event)cudaEventDestroy(event);
      if(completion)cudaEventDestroy(completion);
      throw;
    }
  }
  ~Decoder() { cudaSetDevice(options.gpu); for(auto event:timings) if(event)cudaEventDestroy(event); if(completion)cudaEventDestroy(completion); }
  void mark(int i){if(profile)check_cuda(cudaEventRecord(timings[i],gpu.stream),"record profile event");}
};

std::unique_ptr<Batch> read_batch(int fd, const uint64_t* offsets, const uint64_t* lengths, size_t count, const char* metadata, size_t metadata_size, const int64_t* positions) {
  const char* metadata_end=metadata+metadata_size;
  auto expected_text=[&](){const char* start=metadata;while(metadata<metadata_end && *metadata)++metadata;if(metadata==metadata_end)fail("truncated expected metadata");std::string result(start,metadata-start);++metadata;return result;};
  auto batch = std::make_unique<Batch>();
  if (!count) return batch;
  if (fd<0) fail("invalid BGEN descriptor");
  // Every compressed payload is smaller than its source record. Reserve once
  // (including at most 15 alignment bytes per record) to avoid repeated growth
  // and whole-buffer copies while parsing a multi-megabyte batch.
  uint64_t packed_capacity=0;
  for(size_t i=0;i<count;++i){
    if(lengths[i]>(1ULL<<32) || packed_capacity>(1ULL<<32)-lengths[i])
      fail("invalid/oversized total BGEN batch");
    packed_capacity+=lengths[i];
  }
  if(count>((1ULL<<32)-packed_capacity)/15)fail("BGEN batch alignment exceeds safety bound");
  batch->compressed.reserve(static_cast<size_t>(packed_capacity+15*count));
  batch->records.reserve(count);
  size_t raw_total=0;
  for (size_t first=0;first<count;) {
    size_t end=first+1;
    uint64_t stop=offsets[first]+lengths[first];
    while(end<count && offsets[end]==stop) { stop+=lengths[end]; ++end; }
    if(stop<offsets[first] || stop-offsets[first] > (1ULL<<32)) fail("invalid/oversized BGEN read batch");
    std::vector<uint8_t> block(stop-offsets[first]);
    pread_full(fd,block.data(),block.size(),offsets[first]);
    for(size_t i=first;i<end;++i) {
      const uint8_t* record=block.data()+offsets[i]-offsets[first];
      size_t size=lengths[i], p=0;
      auto skip=[&](size_t n){if(n>size-p) fail("truncated BGEN record"); p+=n;};
      auto text16=[&](){if(size-p<2)fail("truncated text length");size_t n=le16(record+p);skip(2);size_t start=p;skip(n);return std::string(reinterpret_cast<const char*>(record+start),n);};
      std::string variant_id=text16(),rsid=text16(),chromosome=text16();
      if(size-p<4)fail("truncated position");uint32_t position=le32(record+p);skip(4);
      std::string expected_chrom=expected_text(),expected_id=expected_text(),expected_a1=expected_text(),expected_a2=expected_text();
      if(chromosome!=expected_chrom || rsid!=expected_id || int64_t(position)!=positions[i])fail("BGI/BGEN metadata mismatch");
      if(size-p<2) fail("truncated allele count");
      uint16_t alleles=le16(record+p);skip(2);
      if(alleles!=2) fail("direct BGEN decoder requires biallelic variants");
      for(int field=0;field<2;++field){if(size-p<4) fail("truncated allele length");size_t n=le32(record+p);skip(4);size_t start=p;skip(n);std::string allele(reinterpret_cast<const char*>(record+start),n);if(allele!=(field==0?expected_a1:expected_a2))fail("BGI/BGEN allele mismatch");}
      if(size-p<8) fail("truncated compressed header");
      uint32_t compressed_length=le32(record+p);skip(4);
      uint32_t raw_size=le32(record+p);skip(4);
      if(compressed_length!=size-p+4 || size-p<6) fail("BGEN compressed length mismatch");
      uint8_t cmf=record[p], flg=record[p+1];
      if((cmf&15)!=8 || (cmf>>4)>7 || (((uint16_t(cmf)<<8)|flg)%31)!=0 || (flg&32)) fail("unsupported zlib header");
      GpuRecord r;
      r.comp_offset=align_up(batch->compressed.size(),16);
      r.comp_size=size-p-6;
      r.raw_offset=align_up(raw_total,16);r.raw_size=raw_size;
      raw_total=r.raw_offset+r.raw_size;
      if(raw_total>(1ULL<<34)) fail("BGEN inflated batch exceeds 16 GiB safety bound; reduce chunk size");
      r.expected_adler=be32(record+size-4);
      batch->compressed.resize(r.comp_offset,0);
      batch->compressed.insert(batch->compressed.end(),record+p+2,record+size-4);
      batch->records.push_back(r);
    }
    first=end;
  }
  return batch;
}

void decode_batch(Decoder& decoder, const Batch& batch, float* output, int* statuses, cudaStream_t allocation_stream) {
  auto& g=decoder.gpu; const auto& o=decoder.options;
  check_cuda(cudaSetDevice(o.gpu),"select decoder GPU");
  cudaEvent_t allocation_ready;
  check_cuda(cudaEventCreateWithFlags(&allocation_ready,cudaEventDisableTiming),"allocation event");
  check_cuda(cudaEventRecord(allocation_ready,allocation_stream),"allocation stream ready");
  check_cuda(cudaStreamWaitEvent(g.stream,allocation_ready,0),"wait allocator dependencies");
  check_cuda(cudaEventDestroy(allocation_ready),"destroy allocation event");
  size_t n=batch.records.size(), raw_total=0,max_raw=0;
  if(!n)return;
  for(const auto& r:batch.records){raw_total=std::max(raw_total,r.raw_offset+r.raw_size);max_raw=std::max(max_raw,r.raw_size);}
  ensure_device(g.d_comp,g.comp_cap,batch.compressed.size());
  ensure_device(g.d_raw,g.raw_cap,raw_total);
  ensure_pinned(g.h_comp,g.h_comp_cap,batch.compressed.size());
  ensure_pinned(g.h_parse_status,g.h_status_cap,n);
  ensure_record_buffers(g,n);
  std::memcpy(g.h_comp,batch.compressed.data(),batch.compressed.size());
  decoder.mark(0);
  check_cuda(cudaMemcpyAsync(g.d_comp,g.h_comp,batch.compressed.size(),cudaMemcpyHostToDevice,g.stream),"compressed H2D");
  std::vector<const void*> cp(n);std::vector<void*> rp(n);
  std::vector<size_t> cs(n),rs(n);std::vector<uint32_t> ad(n);
  for(size_t i=0;i<n;++i){const auto& r=batch.records[i];cp[i]=g.d_comp+r.comp_offset;rp[i]=g.d_raw+r.raw_offset;cs[i]=r.comp_size;rs[i]=r.raw_size;ad[i]=r.expected_adler;}
  check_cuda(cudaMemcpyAsync(g.d_comp_ptrs,cp.data(),n*sizeof(void*),cudaMemcpyHostToDevice,g.stream),"compressed pointers");
  check_cuda(cudaMemcpyAsync(g.d_raw_ptrs,rp.data(),n*sizeof(void*),cudaMemcpyHostToDevice,g.stream),"raw pointers");
  check_cuda(cudaMemcpyAsync(g.d_comp_sizes,cs.data(),n*sizeof(size_t),cudaMemcpyHostToDevice,g.stream),"compressed sizes");
  check_cuda(cudaMemcpyAsync(g.d_raw_caps,rs.data(),n*sizeof(size_t),cudaMemcpyHostToDevice,g.stream),"raw sizes");
  check_cuda(cudaMemcpyAsync(g.d_adler,ad.data(),n*sizeof(uint32_t),cudaMemcpyHostToDevice,g.stream),"checksum values");
  auto opts=nvcompBatchedDeflateDecompressDefaultOpts;opts.backend=NVCOMP_DECOMPRESS_BACKEND_CUDA;
  size_t temp=0;
  check_nvcomp(nvcompBatchedDeflateDecompressGetTempSizeAsync(n,max_raw,opts,&temp,raw_total),"nvCOMP workspace");
  ensure_device(g.d_temp,g.temp_cap,temp);
  // Device capacities as high-water marks. `d_temp` is nvCOMP's decompression
  // scratch, sized by nvcompBatchedDeflateDecompressGetTempSizeAsync at
  // runtime: it is opaque to any static model, varies with nvCOMP version and
  // GPU, and is therefore the one BGEN memory term that must be measured on
  // the machine rather than predicted. Reporting it is what lets the shipped
  // calculator stay honest about device-decoded formats.
  decoder.totals[10]=std::max(decoder.totals[10],double(g.temp_cap));
  decoder.totals[11]=std::max(decoder.totals[11],double(g.comp_cap));
  decoder.totals[12]=std::max(decoder.totals[12],double(g.raw_cap));
  decoder.totals[13]=std::max(decoder.totals[13],double(temp));
  decoder.totals[14]=std::max(decoder.totals[14],double(n));
  decoder.mark(1);
  check_nvcomp(nvcompBatchedDeflateDecompressAsync(g.d_comp_ptrs,g.d_comp_sizes,g.d_raw_caps,g.d_actual_sizes,n,g.d_temp,temp,g.d_raw_ptrs,opts,reinterpret_cast<nvcompStatus_t*>(g.d_nv_status),g.stream),"nvCOMP DEFLATE");
  decoder.mark(2);
  decode_layout2_kernel<<<n,256,0,g.stream>>>(reinterpret_cast<const uint8_t*const*>(g.d_raw_ptrs),g.d_actual_sizes,g.d_raw_caps,g.d_nv_status,g.d_sample_indices,o.n_samples,o.n_bgen_samples,output,g.d_parse_status,n);
  check_cuda(cudaGetLastError(),"Layout2 decode");
  decoder.mark(3);
  if (o.verify_adler) {
    adler32_kernel<<<n,256,0,g.stream>>>(reinterpret_cast<const uint8_t*const*>(g.d_raw_ptrs),g.d_actual_sizes,g.d_adler,g.d_nv_status,g.d_parse_status,n);
    check_cuda(cudaGetLastError(),"Adler32 validation");
  }
  decoder.mark(4);
  check_cuda(cudaMemcpyAsync(g.h_parse_status,g.d_parse_status,n*sizeof(int),cudaMemcpyDeviceToHost,g.stream),"status D2H");
  decoder.mark(5);
  if(decoder.completion){
    check_cuda(cudaEventRecord(decoder.completion,g.stream),"record completion");
    check_cuda(cudaEventSynchronize(decoder.completion),"blocking decode completion");
  }else check_cuda(cudaStreamSynchronize(g.stream),"decode completion");
  if(decoder.profile){
    decoder.totals[1]+=1;decoder.totals[2]+=batch.compressed.size();
    for(const auto& record:batch.records)decoder.totals[3]+=record.raw_size;
    decoder.totals[4]+=double(n)*o.n_samples*sizeof(float);
    for(int i=0;i<5;++i){float ms=0;check_cuda(cudaEventElapsedTime(&ms,decoder.timings[i],decoder.timings[i+1]),"profile elapsed");decoder.totals[5+i]+=ms;}
  }
  std::memcpy(statuses,g.h_parse_status,n*sizeof(int));
}
__global__ void tg_bgen_probe_kernel(int* out){ if (threadIdx.x == 0) *out = 1; }
} // namespace

extern "C" {
int tg_bgen_abi_version(){return 4;}

/* Can this library's kernels actually launch on `device`?
 *
 * The library loads on any card and ctypes resolves every symbol, but it
 * carries code only for the architectures it was built for -- so an
 * unsupported GPU fails at *launch*, with "no kernel image is available for
 * execution on the device", deep inside a scan. There was no way to ask the
 * question except by attempting a real decode, which needs a real batch, so
 * `resolve_decode_backend` caught only the load-time ImportError and let the
 * launch failure escape.
 *
 * This launches a kernel that does nothing, which is the cheapest question
 * that still exercises the thing that fails. Returns 1 when it ran, 0
 * otherwise; never throws, and clears any error it provoked so the caller's
 * next CUDA call is not poisoned by it.
 */
int tg_bgen_probe(int device){
  int previous = 0;
  if (cudaGetDevice(&previous) != cudaSuccess) return 0;
  if (cudaSetDevice(device) != cudaSuccess) { cudaGetLastError(); return 0; }
  int* slot = nullptr;
  int ok = 0;
  if (cudaMalloc(&slot, sizeof(int)) == cudaSuccess) {
    tg_bgen_probe_kernel<<<1, 32>>>(slot);
    ok = (cudaGetLastError() == cudaSuccess && cudaDeviceSynchronize() == cudaSuccess);
    cudaFree(slot);
  }
  cudaGetLastError();
  cudaSetDevice(previous);
  return ok ? 1 : 0;
}
void tg_bgen_set_verify_adler(void* d,int on){static_cast<Decoder*>(d)->options.verify_adler = on != 0;}
const char* tg_bgen_error(){return last_error.c_str();}
void* tg_bgen_create(int device,size_t n,const int32_t* selection,size_t selected){try{return new Decoder(device,n,selection,selected);}catch(const std::exception& e){last_error=e.what();return nullptr;}}
void tg_bgen_destroy(void* p){delete static_cast<Decoder*>(p);}
void tg_bgen_profile(void* p,double* values){auto* decoder=static_cast<Decoder*>(p);std::copy(decoder->totals,decoder->totals+16,values);}
void* tg_bgen_read(int fd,const uint64_t* offsets,const uint64_t* lengths,size_t count,const char* metadata,size_t metadata_size,const int64_t* positions){try{return read_batch(fd,offsets,lengths,count,metadata,metadata_size,positions).release();}catch(const std::exception& e){last_error=e.what();return nullptr;}}
void tg_bgen_free_batch(void* p){delete static_cast<Batch*>(p);}
int tg_bgen_decode(void* d,void* b,void* out,int* statuses,void* stream){try{decode_batch(*static_cast<Decoder*>(d),*static_cast<Batch*>(b),static_cast<float*>(out),statuses,static_cast<cudaStream_t>(stream));return 0;}catch(const std::exception& e){last_error=e.what();auto* decoder=static_cast<Decoder*>(d);cudaSetDevice(decoder->options.gpu);cudaStreamSynchronize(decoder->gpu.stream);return -1;}}
}
