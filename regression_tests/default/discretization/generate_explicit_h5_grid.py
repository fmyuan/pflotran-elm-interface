
import sys
path_to_src_python = '../../../src/python/unstructured_grid'
sys.path.append(path_to_src_python)
import convert_explicit_ascii_to_h5 as mod

if __name__ == "__main__":
  ret_code = mod.convert_grid_explicit_ascii_to_h5("mixed.uge","mixed_explicit.h5")
