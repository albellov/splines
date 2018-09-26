import matplotlib.pyplot as plt
import sys


def reading(filename):
    vals = [[], []]

    with open(filename, 'r') as f:
        labels = f.readline().split(',')
        labels[-1] = labels[-1].replace('\n', '')

        for line in f:
            spl_line = line.split(',')
            vals[0].append(float(spl_line[0]))
            vals[1].append(float(spl_line[1]))

    return vals, labels


def draw_2d(data_filename, result_filename, result_knots_filename):
    vals, labels = reading(data_filename)
    plt.plot(vals[0], vals[1], 'r.', label='data')

    vals, labels = reading(result_filename)
    plt.plot(vals[0], vals[1], 'b', label='spline')

    vals, labels = reading(result_knots_filename)
    plt.plot(vals[0], vals[1], 'gs', label='knots')

    plt.title('Approximation of spline')
    plt.xlabel('X', fontsize=12)
    plt.ylabel('Y', fontsize=12)

    plt.grid()
    plt.legend()

    plt.show()


def main():
    try:
        dim = int(sys.argv[1])
        data_filename = sys.argv[2]
        result_filename = sys.argv[3]
        result_knots_filename = sys.argv[4]

    except IndexError:
        print("Args: <dimension> <data filename> <result filename> <result knots filename>")
        return
    except ValueError:
        print("Dimension should be int")
        return
    else:
        if dim == 2:
            draw_2d(data_filename, result_filename, result_knots_filename)
        elif dim == 3:
            pass
        else:
            print("Dimension should be 2 or 3")

main()
