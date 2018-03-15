from math import sqrt

N = 7


class Point3(object):
    def __init__(self, x=0, y=0, z=0):
        self.x = x * 1.0
        self.y = y * 1.0
        self.z = z * 1.0

    def __str__(self):
        return ' '.join([str(self.x), str(self.y), str(self.z)])


def dis(a, b):
    return sqrt(pow(a.x-b.x, 2) + pow(a.y-b.y, 2) + pow(a.z-b.z, 2))


def get_nearest(a, input_list):
    d_list = list((input_list.index(i),dis(a,i)) for i in input_list)
    sorted_list = sorted(d_list, key=lambda e:e[1])
    return sorted(list(i[0] for i in sorted_list[:N]))


if __name__ == "__main__":
    Pointlist = []
    for i in range(-1, 2):
        for j in range(-1, 2):
            for k in range(-1, 2):
                Pointlist.append(Point3(i, j, k))

    RNAGE = 8
    delta = 0.5 / RNAGE

    my_dict={}

    for i in range(0, RNAGE):
        for j in range(0, RNAGE):
            for k in range(0, RNAGE):
                new_point = Point3((2*i+1)*delta, (2*j+1)*delta, (2*k+1)*delta)
                key = get_nearest(new_point, Pointlist)
                key_string = ' '.join(list(str(a) for a in key))
                if key_string not in my_dict:
                    my_dict[key_string] = [str(new_point)]
                else:
                    my_dict[key_string].append(str(new_point))

    print len(my_dict.keys())
    print my_dict.keys()
