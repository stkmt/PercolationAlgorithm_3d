#!/usr/bin/env python

import numpy as np
import sys
import time

#Date 2021/03
#Author Shuta Takimoto

# the lines commented out is gonna be useful if you wanna get more information.

argvs = sys.argv
start_time = time.time()

if len(argvs) != 2:
    print("\n"
          "         PROGRAM DESCRIPTION\n"
          "\n"
          "         This program is going to estimate migration barrier by using a percolation algorithm.\n"
          "         'BVSx_xyzp.dat' (3D model) is necessary to run this program because this program (for calculating\n"
          "         migration barrier) requires it.\n"
          "         COMMAND --> 'python Percolation.py BVSx_xyzp.dat'\n"
          "\n"
          "         A cluster-invade type site percolation algorithm is adopted here.\n"
          "         So it's clearly different from a ravel type site percolation or any other bond percolation algorithms.\n"
          "\n"
          "         The migration barrier in the path calculated by this program will be the difference between the maximum\n"
          "         and the minimum in the path.\n"
          "         Meanwhile, the migration barrier in the lattice will be the difference between the maximum in the path\n"
          "         and the minimum in the lattice.\n"
          "         As you know, 'Em lattice  >  Em path' must be basically seen.\n"
          "\n"
          "         Hope your success. Bye:)\n"
          )

    sys.exit()


class Percolation:

    def __init__(self, direction, xmax, ymax, zmax, potential_map_original):
        self.dire = direction
        self.xmax = xmax
        self.ymax = ymax
        self.zmax = zmax
        self.p_map = np.empty((xmax, ymax, zmax))
        self.add_potential(potential_map_original)


    def add_potential(self, potential_map_original): # initialize potential_map
        for n, p in enumerate(np.ravel(potential_map_original)):
            np.put(self.p_map, [n], p)


    def first_surface(self, potential_map_original): # explore which site is active on the first surface
        surface_list = []
        pote_map = np.empty((self.xmax, self.ymax, self.zmax))
        for n, p in enumerate(np.ravel(potential_map_original)):
            np.put(pote_map, [n], p)

        if self.dire == 'X':
            for y in range(self.ymax):
                for z in range(self.zmax):
                    if self.p_map[0, y, z] != 0:
                        surface_list.append([0, y, z])

        elif self.dire == 'Y':
            for x in range(self.xmax):
                for z in range(self.zmax):
                    if self.p_map[x, 0, z] != 0:
                        surface_list.append([x, 0, z])

        elif self.dire == 'Z':
            for x in range(self.xmax):
                for y in range(self.ymax):
                    if self.p_map[x, y, 0] != 0:
                        surface_list.append([x, y, 0])

        self.p_map_dict = {str(s[0])+','+str(s[1])+','+str(s[2]): pote_map for s in surface_list}

        return surface_list


    def move_next(self, current_list, surface_list, path_list): # move to next site
        next_list = []
        surface_list_out = []

        for candidate, start in zip(current_list, surface_list):
            self.p_map_dict[str(start[0])+','+str(start[1])+','+str(start[2])][candidate[0], candidate[1], candidate[2]] = 0

            if self.dire == 'X':
                if candidate[0] != self.xmax-1:
                    next_site = self.check_next(candidate, start)
                    if not next_site:
                        continue

                    else:
                        next_list.append(next_site)
                        surface_list_out.append(start)

                else:
                    path_list.append('{0},{1},{2}-{3},{4},{5}'.format(start[0], start[1], start[2], candidate[0], candidate[1], candidate[2]))

            if self.dire == 'Y':
                if candidate[1] != self.ymax-1:
                    next_site = self.check_next(candidate, start)
                    if not next_site:
                        continue

                    else:
                        next_list.append(next_site)
                        surface_list_out.append(start)

                else:
                    path_list.append('{0},{1},{2}-{3},{4},{5}'.format(start[0], start[1], start[2], candidate[0], candidate[1], candidate[2]))

            if self.dire == 'Z':
                if candidate[2] != self.zmax-1:
                    next_site = self.check_next(candidate, start)
                    if not next_site:
                        continue

                    else:
                        next_list.append(next_site)
                        surface_list_out.append(start)

                else:
                    path_list.append('{0},{1},{2}-{3},{4},{5}'.format(start[0], start[1], start[2], candidate[0], candidate[1], candidate[2]))

        return next_list, surface_list_out, path_list


    def check_next(self, current_site, start_site): # detect next site with considering periodic boundary condition
        check_box = {}
        next_sites = [[current_site[0]+1, current_site[1], current_site[2]],
                      [current_site[0]-1, current_site[1], current_site[2]],
                      [current_site[0], current_site[1]+1, current_site[2]],
                      [current_site[0], current_site[1]-1, current_site[2]],
                      [current_site[0], current_site[1], current_site[2]+1],
                      [current_site[0], current_site[1], current_site[2]-1]]

        if next_sites[0][0] > self.xmax-1:
            next_sites[0][0] = next_sites[0][0] - (self.xmax-1)
        if next_sites[1][0] < 0:
            next_sites[1][0] = next_sites[1][0] + self.xmax-1
        if next_sites[2][1] > self.ymax-1:
            next_sites[2][1] = next_sites[2][1] - (self.ymax-1)
        if next_sites[3][1] < 0:
            next_sites[3][1] = next_sites[3][1] + self.ymax-1
        if next_sites[4][2] > self.zmax-1:
            next_sites[4][2] = next_sites[4][2] - (self.zmax-1)
        if next_sites[5][2] < 0:
            next_sites[5][2] = next_sites[5][2] + self.zmax-1

        if self.dire == 'X':
            if current_site[0] != 0:
                behind_site = next_sites[1]

            else:
                 behind_site = None

            del next_sites[1]
            for n in next_sites:
                check_box[str(n[0])+','+str(n[1])+','+str(n[2])] = self.p_map_dict[str(start_site[0])+','+str(start_site[1])+','+str(start_site[2])][n[0], n[1], n[2]]

        elif self.dire == 'Y':
            if current_site[1] != 0:
                behind_site = next_sites[3]

            else:
                 behind_site = None

            del next_sites[3]
            for n in next_sites:
                check_box[str(n[0])+','+str(n[1])+','+str(n[2])] = self.p_map_dict[str(start_site[0])+','+str(start_site[1])+','+str(start_site[2])][n[0], n[1], n[2]]

        elif self.dire == 'Z':
            if current_site[2] != 0:
                behind_site = next_sites[5]

            else:
                 behind_site = None

            del next_sites[5]
            for n in next_sites:
                check_box[str(n[0])+','+str(n[1])+','+str(n[2])] = self.p_map_dict[str(start_site[0])+','+str(start_site[1])+','+str(start_site[2])][n[0], n[1], n[2]]

        if list(check_box.values()) == [0, 0, 0, 0, 0]:
            if behind_site != None and self.p_map_dict[str(start_site[0])+','+str(start_site[1])+','+str(start_site[2])][behind_site[0], behind_site[1], behind_site[2]] == 0:
                next_site = None

            else:
                next_site = behind_site

        else:
            min_check_box = self.pick_out_of_dict(check_box)
            next_site = [int(min_check_box.split(',')[0]), int(min_check_box.split(',')[1]), int(min_check_box.split(',')[2])]

        return next_site


    def stock_potential(self, current_list, max_potential_stock_list, min_potential_stock_list,  surface_list, potential_dict): # stock min and max potential in the path
        max_potential_stock_list_out = []
        min_potential_stock_list_out = []

        potential_dict_out = potential_dict

        for maxp, minp, c, s in zip(max_potential_stock_list, min_potential_stock_list, current_list, surface_list):
            if maxp < self.p_map[c[0], c[1], c[2]] and self.p_map[c[0], c[1], c[2]] != 0:
                maxp = self.p_map[c[0], c[1], c[2]]

            if minp > self.p_map[c[0], c[1], c[2]]:
                minp = self.p_map[c[0], c[1], c[2]]

            max_potential_stock_list_out.append(maxp)
            min_potential_stock_list_out.append(minp)
            potential_dict_out[str(s[0])+','+str(s[1])+','+str(s[2])] = [minp, maxp]

        return max_potential_stock_list_out, min_potential_stock_list_out, potential_dict_out


    def pick_out_of_dict(self, dictionary, mode = 'min'): # get a key of the minimum value of a dictionary
        if mode == 'min':
            objective = min(dictionary, key=dictionary.get)
        elif mode == 'max':
            objective = max(dictionary, key=dictionary.get)

        return objective


    def data(self, path_list, potential_dict): # estimate Ea etc..
        diff_potential_dict = {key: value[1]-value[0] for key, value in potential_dict.items()}
        min_potential_dict = {key: value[0] for key, value in potential_dict.items()}
        max_potential_dict = {key: value[1] for key, value in potential_dict.items()}

        if not potential_dict:
            diff_potential_dict = {'1000,1000,1000': 1000}
            min_potential_dict = {'1000,1000,1000': 0}
            max_potential_dict = {'1000,1000,1000': 1000}

        unstable_key = self.pick_out_of_dict(max_potential_dict)
        unstable_value = max_potential_dict[unstable_key]
#        stable_key = self.pick_out_of_dict(min_potential_dict)
#        stable_value = min_potential_dict[stable_key]

#        stable_from_unstable = min_potential_dict[unstable_key]
#        unstable_from_stable = max_potential_dict[stable_key]

        final_start_list = [i.split('-')[0].split(',') for i in path_list]
        final_start_list = [[int(s) for s in i] for i in final_start_list]
        final_goal_list = [i.split('-')[1].split(',') for i in path_list]
        final_goal_list = [[int(s) for s in i] for i in final_goal_list]

        if not path_list:
            final_start_list = [[1000, 1000, 1000]]
            final_goal_list = [[1000, 1000, 1000]]

        min_potential_lattice = np.min(self.p_map)

        diff_key = self.pick_out_of_dict(diff_potential_dict)
        diff_value = diff_potential_dict[diff_key]
#        max_min = unstable_value - stable_from_unstable

#        self.start = [int(unstable_key.split(',')[0]),int(unstable_key.split(',')[1]),int(unstable_key.split(',')[2])]
        self.start = [int(diff_key.split(',')[0]), int(diff_key.split(',')[1]), int(diff_key.split(',')[2])]
        self.goal = final_goal_list[final_start_list.index(self.start)]

        self.diff_ea = diff_value # the smallest difference between high and low
#        self.max_min_ea = max_min # the lowest high, and its low
        self.lattice_ea = unstable_value - min_potential_lattice # the lowest high and the lowest lattice potential

        if self.dire == 'X':
            self.path_ratio = len(path_list) / (self.ymax * self.zmax)

        elif self.dire == 'Y':
            self.path_ratio = len(path_list) / (self.xmax * self.zmax)

        elif self.dire == 'Z':
            self.path_ratio = len(path_list) / (self.xmax * self.ymax)

        print("\n")
        print("Direction: {}\n".format(self.dire))
        print("Start site: {}\n".format(self.start))
        print("Goal site: {}\n".format(self.goal))
        print("Minimum potential in the path: {} eV\n".format(min_potential_dict[diff_key]))
        print("Maximum potential in the path: {} eV\n".format(max_potential_dict[diff_key]))
#        print("Minimum potential: {} eV\n".format(stable_from_unstable))
#        print("Maximum potential: {} eV\n".format(unstable_value))
        print("Migration barrier: {:.5f} eV\n".format(self.diff_ea))
        print("the num of path / surface area: {:.3f}\n".format(self.path_ratio))


    def memory(self, potential_map_original): # memorize the path and potential values in it
        self.path_memo = np.zeros((self.xmax, self.ymax, self.zmax))
        self.p_map = np.empty((self.xmax, self.ymax, self.zmax))
        self.add_potential(potential_map_original)
        self.p_map_dict = {str(self.start[0])+','+str(self.start[1])+','+str(self.start[2]): self.p_map}
        n = 0

        if self.start != [1000, 1000, 1000]:
            while True:
                if n == 0:
                    self.path_memo[self.start[0], self.start[1], self.start[2]] = self.p_map_dict[str(self.start[0])+','+str(self.start[1])+','+str(self.start[2])][self.start[0], self.start[1], self.start[2]]
                    self.p_map_dict[str(self.start[0])+','+str(self.start[1])+','+str(self.start[2])][self.start[0], self.start[1], self.start[2]] = 0
                    next_site = self.check_next(self.start, self.start)
                    current_site = next_site

                else:
                    self.path_memo[current_site[0], current_site[1], current_site[2]] = self.p_map_dict[str(self.start[0])+','+str(self.start[1])+','+str(self.start[2])][current_site[0], current_site[1], current_site[2]]
                    self.p_map_dict[str(self.start[0])+','+str(self.start[1])+','+str(self.start[2])][current_site[0], current_site[1], current_site[2]] = 0
                    next_site = self.check_next(current_site, self.start)
                    current_site = next_site

                if current_site == None:
                    break

                if self.dire == 'X' and current_site[0] == self.xmax-1:
                    break

                elif self.dire == 'Y' and current_site[1] == self.ymax-1:
                    break

                elif self.dire == 'Z' and current_site[2] == self.zmax-1:
                    break

                n = n + 1

        else:
            for i in range(self.xmax*self.ymax*self.zmax):
                np.put(self.path_memo, [i], i)


if __name__=='__main__':

    inpf = argvs[1]

    with open(inpf) as f:
        line = f.readline()
        l = 0
        for i in f:
            if l == 0:
                xmax = int(i.split()[1])
                ymax = int(i.split()[2])
                zmax = int(i.split()[3])

                potential_map = np.empty((xmax, ymax, zmax))

            elif l >= 1:
                np.put(potential_map, [l-1], float(i.split()[3]))

            l = l + 1


    xpath = Percolation('X', xmax, ymax, zmax, potential_map) # ------------------- X direction

    start_list = xpath.first_surface(potential_map)
    path_list = []
    potential_dict = {}
    n = 0

    while True:
        if n == 0:
            next_list, start_list, path_list = xpath.move_next(start_list, start_list, path_list)
            current_list = next_list

            max_potential_list = [-100 for i in range(len(current_list))]
            min_potential_list = [100 for i in range(len(current_list))]
            max_potential_list, min_potential_list, potential_dict = xpath.stock_potential(current_list, max_potential_list, min_potential_list, start_list, potential_dict)

        else:
            next_list, start_list, path_list = xpath.move_next(current_list, start_list, path_list)
            current_list = next_list

            max_potential_list, min_potential_list, potential_dict = xpath.stock_potential(current_list, max_potential_list, min_potential_list, start_list, potential_dict)

        if not current_list:
            break

        n = n + 1

    potential_dict = {key: value for key, value in potential_dict.items() for i in path_list if key == i.split('-')[0]}

    xpath.data(path_list, potential_dict)
    xpath.memory(potential_map)

    ypath = Percolation('Y', xmax, ymax, zmax, potential_map) # ------------------- Y direction

    start_list = ypath.first_surface(potential_map)
    path_list = []
    potential_dict = {}
    n = 0

    while True:
        if n == 0:
            next_list, start_list, path_list = ypath.move_next(start_list, start_list, path_list)
            current_list = next_list

            max_potential_list = [-100 for i in range(len(current_list))]
            min_potential_list = [100 for i in range(len(current_list))]
            max_potential_list, min_potential_list, potential_dict = ypath.stock_potential(current_list, max_potential_list, min_potential_list, start_list, potential_dict)

        else:
            next_list, start_list, path_list = ypath.move_next(current_list, start_list, path_list)
            current_list = next_list

            max_potential_list, min_potential_list, potential_dict = ypath.stock_potential(current_list, max_potential_list, min_potential_list, start_list, potential_dict)

        if not current_list:
            break

        n = n + 1

    potential_dict = {key: value for key, value in potential_dict.items() for i in path_list if key == i.split('-')[0]}

    ypath.data(path_list, potential_dict)
    ypath.memory(potential_map)

    zpath = Percolation('Z', xmax, ymax, zmax, potential_map) # ------------------- Z direction

    start_list = zpath.first_surface(potential_map)
    path_list = []
    potential_dict = {}
    n = 0

    while True:
        if n == 0:
            next_list, start_list, path_list = zpath.move_next(start_list, start_list, path_list)
            current_list = next_list

            max_potential_list = [-100 for i in range(len(current_list))]
            min_potential_list = [100 for i in range(len(current_list))]
            max_potential_list, min_potential_list, potential_dict = zpath.stock_potential(current_list, max_potential_list, min_potential_list, start_list, potential_dict)

        else:
            next_list, start_list, path_list = zpath.move_next(current_list, start_list, path_list)
            current_list = next_list

            max_potential_list, min_potential_list, potential_dict = zpath.stock_potential(current_list, max_potential_list, min_potential_list, start_list, potential_dict)

        if not current_list:
            break

        n = n + 1

    potential_dict = {key: value for key, value in potential_dict.items() for i in path_list if key == i.split('-')[0]}

    zpath.data(path_list, potential_dict)
    zpath.memory(potential_map)

    with open('out.perco', 'w') as f:
        f.write("# minPot: {} eV\n".format(np.min(potential_map)))
        f.write("# Mesh: X: {0}, Y: {1}, Z: {2}\n".format(xmax, ymax, zmax))
        f.write("# Em(eV vs. in latiice) Em(eV vs. in path) path_ratio direction\n".format(np.min(potential_map)))
        f.write("  {0:<21.5f} {1:<18.5f} {2:<10.3f} {3:<9}\n".format(xpath.lattice_ea, xpath.diff_ea, xpath.path_ratio, xpath.dire))
        f.write("  {0:<21.5f} {1:<18.5f} {2:<10.3f} {3:<9}\n".format(ypath.lattice_ea, ypath.diff_ea, ypath.path_ratio, ypath.dire))
        f.write("  {0:<21.5f} {1:<18.5f} {2:<10.3f} {3:<9}\n".format(zpath.lattice_ea, zpath.diff_ea, zpath.path_ratio, zpath.dire))

    with open('out.xpath', 'w') as f:
        for i in xpath.path_memo.ravel():
            f.write(str(i) + "\n")

    with open('out.ypath', 'w') as f:
        for i in ypath.path_memo.ravel():
            f.write(str(i) + "\n")

    with open('out.zpath', 'w') as f:
        for i in zpath.path_memo.ravel():
            f.write(str(i) + "\n")

    print("\n")
    print("---------------------")
    print("output: out.perco")
    print("output: out.xpath")
    print("output: out.ypath")
    print("output: out.zpath")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print("\n")
    print("Elapsed time: {} s\n".format(int(elapsed_time)))
