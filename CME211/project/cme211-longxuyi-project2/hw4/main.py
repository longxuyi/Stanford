import sys
import truss


if len(sys.argv) < 3:
    print('Usage:')
    print('  python3 main.py <joints_file> <beams_file> <optional_plot_file>')
    sys.exit(0)

joints_file = sys.argv[1]
beams_file = sys.argv[2]

#if optional input is not given
if len(sys.argv) == 3:
    plot_file = ""

#if optional input is given
if len(sys.argv) == 4:
    plot_file = sys.argv[3]

a = truss.Truss(joints_file, beams_file, plot_file)
print(a)