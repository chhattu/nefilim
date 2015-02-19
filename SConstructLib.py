#!/usr/bin/python

import os

# Compilation_types = enum(CPU=1, MPI=2)

def get_rosetta_shared_lib_path(rosetta_home, type_of_compilation):
  rosetta_build_path      = 'build/src/release/' # join does not take suffix startswith /
  rosetta_shared_lib_path = os.path.join(rosetta_home, rosetta_build_path)
  if os.path.isdir(rosetta_shared_lib_path):
    for root, dirs, files in os.walk(rosetta_shared_lib_path):
      for libfile in ( file for file in files if file.startswith("lib")): 
        if type_of_compilation.lower() == "mpi":
          if root.count("mpi")  > 0: rosetta_shared_lib_path = root 
        elif type_of_compilation.lower() == "single":
          if root.count("mpi") <= 0: rosetta_shared_lib_path = root 
    # print "Debug: shared library path " + rosetta_shared_lib_path
    return rosetta_shared_lib_path 
  else:
    return ""
  
def get_platform_path(rosetta_home, platform_type):
  if len(rosetta_home.strip()) <= 0: return ""
  if not os.path.isdir(rosetta_home): return ""
  if len(platform_type.strip()) <= 0: return ""
  platform_type = platform_type.lower()
  if platform_type.strip() == "darwin": platform_type = "macos" 
  intermediate_path = os.path.join("src/platform/" + platform_type)
  fullpath = os.path.join(rosetta_home, intermediate_path)
  # print "Debug: platform path " + fullpath  # DEBUG
  if os.path.isdir(fullpath): return "/" + intermediate_path
  else: return ""
  

