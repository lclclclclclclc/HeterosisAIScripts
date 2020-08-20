# -*- coding: utf-8 -*-
import random
import numpy as np

offspring_sexes = ['F'] * 50 + ['M'] * 50
random.shuffle(offspring_sexes)

parent_subpops = [3] * 90 + [2] * 10
random.shuffle(parent_subpops)

children = []
while len(offspring_sexes) > 0:
    child_sex = offspring_sexes.pop()

    accept_child = False

    while accept_child is False:
        parent_subpop = parent_subpops.pop()
        if parent_subpop == 3:
            accept_child = True
        elif child_sex == 'M':
            accept_child = True
        parent_subpops.append(parent_subpop)
        random.shuffle(parent_subpops)

    children.append((child_sex, parent_subpop))


migrant_females = 0
migrant_males = 0
local_females = 0
local_males = 0

for c in children:
    if c == ('F', 3):
        local_females += 1
    elif c == ('F', 2):
        migrant_females += 1
    elif c == ('M', 3):
        local_males += 1
    elif c == ('M', 2):
        migrant_males += 1

print(f"Achieved {migrant_females}:{migrant_males} ({migrant_females / migrant_males}) sex-bias ratio,")
print(f"\twith {migrant_males + local_males}% overall male ratio.")
print(f"Maintained {migrant_females + migrant_males}% introgression ratio.")
