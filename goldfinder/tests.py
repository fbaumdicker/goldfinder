import os
import filecmp
import shutil

# Inspired from: https://janakiev.com/blog/python-filecmp/
def gather_files(path): 
    files = []
    if not path.endswith('/'):
        path += '/'
    for (dirpath, dirnames, filenames) in os.walk(path):
        for f in filenames:
            files.append(os.path.join(dirpath, f)[len(path):])
    return files

def perform_tests():

    # cd into goldfinder's dir
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    os.chdir('..')

    # create tmp folder
    path_tmp = './tmp/'
    if not os.path.exists(path_tmp):
        os.makedirs(path_tmp)
    else:
        exit('Please delete the folder "goldfinder/tmp/" manually before running tests again')


    # Run tests
    total_num_tests = 3
    print(f'(1/{total_num_tests}) Running tab test...')
    os.system("python goldfinder/goldfinder.py -t example_files/example.nwk -i example_files/tab_mini_"
            "example.tsv --seed 123 -k example_files/known_assocs.tsv -o tmp/tab -m example_file"
            f"s/metadata.csv -f tab >tmp/tab.log 2>&1")

    print(f'(2/{total_num_tests}) Running roary test...')
    os.system("python goldfinder/goldfinder.py -t example_files/example.nwk -i example_files/roary_min"
            "i_example.csv --seed 123 -k example_files/known_assocs.tsv -o tmp/roary -m example_fi"
            f"les/metadata.csv -f roary >tmp/roary.log 2>&1")

    print(f'(3/{total_num_tests}) Running matrix test...')
    os.system("python goldfinder/goldfinder.py -t example_files/example.nwk -i example_files/matrix_mi"
            "ni_example.csv --seed 123 -k example_files/known_assocs.tsv -o tmp/matrix -m example_fi"
            f"les/metadata.csv >tmp/matrix.log 2>&1")

    # compare directories
    path_exp = './example_files/expected_test_results/'

    match, mismatch, error = filecmp.cmpfiles(path_tmp, path_exp, common = gather_files(path_exp), shallow = False)

    if error:
        print('Tests failed. Following files could not be compared:')
        for file in error:
            print(file)
        print('Please check files and then delete golfinder/tmp folder manually.')

    if mismatch:
        print('Tests failed. Following files yielded different results:')
        for file in mismatch:
            print(file)
        print('Please check files and then delete golfinder/tmp folder manually.')

    if not error and not mismatch:
        print('Tests successful!')
        shutil.rmtree(path_tmp)
