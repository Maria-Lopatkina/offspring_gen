import random

diiiict = {"a": 8, "b": 9, "c": 10}
s = 8 + 9 + 10
ss = 0
while ss != s:
    rand = random.choice(list(diiiict.keys()))
    if diiiict[rand] != 0:
        diiiict[rand] -= 1
        ss += 1
print(diiiict)
