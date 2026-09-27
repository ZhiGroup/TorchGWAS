"""Match measured environment-ready and historical process timing boundaries."""
import copy

def timing_profiles(profiles,boundary='environment-ready'):
 if boundary not in ['environment-ready','process']:raise ValueError('Unknown timing boundary')
 result=copy.deepcopy(profiles)
 for name,profile in result.items():
  profile['timing_boundary']=boundary if name=='torchGWAS' else 'process'
 return result
