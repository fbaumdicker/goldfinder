

#Terminal score (1)
def terminal_score(p,g):
    score = p * g + (1-p)*(1-g) - (1-p)*g - p*(1-g)
    print(score)

terminal_score(1,0)
terminal_score(0,1)
terminal_score(1,1)
terminal_score(0,0)
print("-----------------")
#Simultaneous score (2)
def simu_score(panc,pdesc,ganc,gdes):
    score = (panc - pdesc)*(ganc - gdes)
    return score

for x in range (2):
    for y in range(2):
        for z in range(2):
            for w in range(2):
                score = simu_score(x,y,z,w)
                print(x,y,z,w,"score:", score)

print("-----------------")

#Subsequent score (3)
def sub_score (panc, pdes,ganc,gdes):
    score = 4/3*panc*ganc + 2/3*panc*gdes + 2/3*pdes*ganc + 4/3*pdes*gdes - panc-pdes-ganc-gdes +1
    return score

for x in range (2):
    for y in range(2):
        for z in range(2):
            for w in range(2):
                score = sub_score(x,y,z,w)
                print(x,y,z,w,"score:", score)