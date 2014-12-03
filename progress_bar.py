import sys
import time

toolbar_width = 75

sys.stdout.write("[%s]" % (" " * toolbar_width))
sys.stdout.flush()
sys.stdout.write("\b" * (toolbar_width+1))

for i in xrange(100):
#   time.sleep(0.1)
    num = 100/10
    if i % num == 0:
        sys.stdout.write("-")
        sys.stdout.flush()
sys.stdout.write("\n")    
