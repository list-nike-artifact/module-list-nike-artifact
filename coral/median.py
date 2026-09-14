# Reduces the output of CORAL's own `test_nike_*` benchmark to medians.
#
#     test_nike_p_500 --single --iterations=1000 | tr '\r' '\n' | python3 median.py [GHz]
import re, sys, statistics
ghz = float(sys.argv[1]) if len(sys.argv) > 1 else 4.05
keygen, derive = [], []
for line in sys.stdin:
    m = re.findall(r"This:\s+([\d.]+) ms", line)
    if len(m) == 2:
        keygen.append(float(m[0]))
        derive.append(float(m[1]))
if not keygen:
    sys.exit("no timing lines found on stdin")
print("runs %d  keygen median %.3f ms (%.2f Mcyc)  derive median %.3f ms (%.2f Mcyc)"
      % (len(keygen), statistics.median(keygen), statistics.median(keygen) * ghz,
         statistics.median(derive), statistics.median(derive) * ghz))
