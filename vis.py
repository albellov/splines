import matplotlib.pyplot as plt
import sys


def reading(filename):
    vals = [[],[]]

    with open(filename,'r') as f:

        labels = f.readline().split(',')
        labels[-1] = labels[-1].replace('\n', '')

        for line in f:
            spl_line = line.split(',')
            vals[0].append(float(spl_line[0]))
            vals[1].append(float(spl_line[1]))

    return vals, labels

def main():
    try:
        data_filename = sys.argv[1]
        result_filename = sys.argv[2]
    except:
        print ("Args: <data filename> <in filename>")
        raise SystemExit(0)

    vals, labels = reading(data_filename)
    plt.plot(vals[0], vals[1], 'r.', label='data')

    vals, labels = reading(result_filename)
    plt.plot(vals[0], vals[1], 'b', label='spline')

    plt.title('Approximation of spline')
    plt.xlabel('X', fontsize = 12)
    plt.ylabel('Y', fontsize = 12)
    plt.grid()
    plt.legend()

    plt.show()

main()