import scripts.sys
tom_in = scripts.sys.argv[1]
tom_out = scripts.sys.argv[2]
a = []
ranked = []


def read_files(a):
    with open(tom_in) as txt:
        for line in txt:
            if line.startswith("#") or line == "\n":
                continue
            else:
                c = map(int, line.split()[1:])
                c.insert(0, line.split()[0])
                a += [c]
    return a


def add(an, s):
    for i in range(len(an)):
        if len(an[i]) > 9:
            an[i].pop()
        an[i].insert(9, sum(an[i][s:s+3]))
    an.sort(key=lambda x: x[9], reverse=True)
    return an


def tmp_list(a):
    my_list = []
    for i in a:
        my_list.append(i[9])
    return my_list


def rank(an):
    ranked.append(an[0][0])
    an.pop(0)
    #print a[0]


def isempty(new):
    if len([i for i, x in enumerate(new) if new.count(x) > 1]):
        return True
    else:
        return False


def is_dup(new):
    if len([i for i, x in enumerate(new) if new.count(x) > 1])>0:
        return True
    else:
        return False


def test(a):
    rem = tmp_list(a).count(tmp_list(a)[0])
    if rem == 1:
        rank(a)
    else:
        new = a[:rem]
        pop(rem, a)
        return new


def pop(rem, a):
    for i in range(rem):
        a.pop()


def cond(s, h, rem, new):
    if s > 1:
        #print s
        new = add(new, s-1)
        s -= 1
    else:
        if h <= 7:
            #print h
            new = add(new, h+1)
            h += 1
    return new, s, h

a = read_files(a)
s = 3
add(a, s)
flag = 1
#s,h=repeat(a,s,h)
#print a
while len(a) != 0:
    #print a
    rem = tmp_list(a).count(tmp_list(a)[0])
    #print rem
    if rem == 1:
        rank(a)
    else:
        new = a[:rem]
        for i in range(rem):
            a.pop(0)
        while len(new) != 0:
            new = add(new, s-1)
            rem = tmp_list(new).count(tmp_list(new)[0])
            if rem == 1:
                rank(new)
            else:
                new = add(new, s-2)
                rem = tmp_list(new).count(tmp_list(new)[0])
                if rem == 1:
                    rank(new)
                else:
                    new = add(new, s+1)
                    rem = tmp_list(new).count(tmp_list(new)[0])
                    if rem == 1:
                        rank(new)
                    else:
                        new = add(new, s+2)
                        #print "here"
                        rem = tmp_list(new).count(tmp_list(new)[0])
                        if rem == 1:
                            rank(new)
                            rem = tmp_list(new).count(tmp_list(new)[0])
                        else:
                            new = add(new, s+3)
                            rem = tmp_list(new).count(tmp_list(new)[0])
                            if rem == 1:
                                rank(new)
                            else:
                                new = add(new, s+4)
                                rem = tmp_list(new).count(tmp_list(new)[0])
                                if rem == 1:
                                    rank(new)
                                else:
                                    for i in new:
                                        rank(new)
with open(tom_out, "w") as out:
    for i in ranked:
        out.write(i+"\n")