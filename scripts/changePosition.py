import sys
import glob

newx, newy, newz = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])

filelist = glob.glob('dist_density_*.txt')

# Read first file
with open(filelist[0]) as fh:
    data = fh.readlines()

data = [line.strip().split() for line in data]
idens = int(data[0][0])
x0 = float(data[1][0])
y0 = float(data[2][0])
z0 = float(data[3][0])

lines = []
lines.append('{:<8d} {}\n'.format(idens, 'idens'))
lines.append('{:<8.3f} {}\n'.format(newx, 'x0'))
lines.append('{:<8.3f} {}\n'.format(newy, 'y0'))
lines.append('{:<8.3f} {}\n'.format(newz, 'z0'))

for fn in filelist:
    with open(fn, 'w') as fh:
        fh.writelines(lines)
