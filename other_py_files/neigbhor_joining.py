from skbio import alignment, DNA
from Bio.Phylo.TreeConstruction import _Matrix, DistanceCalculator, DistanceMatrix, DistanceTreeConstructor
from Bio import AlignIO



names = ['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9']
matrix = [[0],
         [0.44084, 0],
         [1.64791843, 0.44084, 0],
         [0.44084, 0, 0.44084, 0],
         [0.44084, 0, 0.44084, 0, 0],
         [1.64791843, 0.44084, 1.64791843, 0.44084, 0.44084, 0],
         [0.44084, 0, 0.44084, 0, 0, 0.44084, 0],
         [1.64791843, 0.44084, 0, 0.44084, 0.44084, 1.64791843, 0.44084, 0],
         [1.64791843, 0.44084, 0, 0.44084, 0.44084, 1.64791843, 0.44084, 0, 0]]
dm = DistanceMatrix(names, matrix)


constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)
print(tree)
"""
Tree(rooted=False)
    Clade(branch_length=0, name='Inner3')
        Clade(branch_length=-0.07943417000000008, name='s4')
        Clade(branch_length=0.9033933900000001, name='s2')
        Clade(branch_length=0.22233543249999993, name='Inner2')
            Clade(branch_length=0.6016237824999999, name='s5')
            Clade(branch_length=0.22233543249999999, name='Inner1')
                Clade(branch_length=0.5275119716666667, name='s3')
                Clade(branch_length=-0.08667197166666674, name='s1')
"""

"""
Tree(rooted=False)
    Clade(branch_length=0, name='Inner7')
        Clade(branch_length=0.0, name='Inner6')
            Clade(branch_length=0.09577980374999998, name='Inner3')
                Clade(branch_length=0.594087686, name='s1')
                Clade(branch_length=-0.15324768599999994, name='s2')
            Clade(branch_length=0.09577980374999998, name='Inner4')
                Clade(branch_length=-0.1436697056250001, name='s4')
                Clade(branch_length=0.5845097056250002, name='Inner2')
                    Clade(branch_length=0.0, name='s9')
                    Clade(branch_length=0.0, name='Inner1')
                        Clade(branch_length=0.0, name='s8')
                        Clade(branch_length=0.0, name='s3')
        Clade(branch_length=0.09577980374999998, name='Inner5')
            Clade(branch_length=0.568546405, name='s6')
            Clade(branch_length=-0.12770640499999997, name='s5')
        Clade(branch_length=-0.09577980374999998, name='s7')
"""