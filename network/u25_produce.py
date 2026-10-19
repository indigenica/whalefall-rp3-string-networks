"""U25 (docs/preregistros/ESPEC_U25.md): the dtau = 0.125 networks alone.  The production queue lives in u24_produce.py,
which interleaves them after each U24 seed (so every pair exists early); this runs only the U25 entries (skipping those done).
    python -B network/u25_produce.py"""
import sys
import u24_produce

if __name__ == "__main__":
    u24_produce.main(sys.argv + ["--only", "U25"])
