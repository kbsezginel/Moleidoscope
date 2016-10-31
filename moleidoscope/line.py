from moleidoscope.linker import Linker


class Line:
    """ Line class."""
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        self.name = str(p1) + '_' + str(p2)
        self.vec = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]

    def grid(self, size):
        line_coors = []
        x, y, z = self.p1
        for s in range(size):
            x += self.vec[0] / size
            y += self.vec[1] / size
            z += self.vec[2] / size
            line_coors.append([x, y, z])
        return line_coors

    def to_linker(self, atom_type='F', size=10):
        line_coors = self.grid(size)
        line_names = [atom_type] * len(line_coors)
        l_line = Linker()
        l_line.name = self.name
        l_line.atom_coors = line_coors
        l_line.atom_names = line_names
        return l_line

    def xyz(length=30, size=10):
        x_line = Line([0, 0, 0], [length, 0, 0])
        y_line = Line([0, 0, 0], [0, length, 0])
        z_line = Line([0, 0, 0], [0, 0, length])

        lx = x_line.to_linker(atom_type='O', size=size)
        ly = y_line.to_linker(atom_type='N', size=size)
        lz = z_line.to_linker(atom_type='C', size=size)
        lxyz = lx.join(ly, lz)
        return lxyz
