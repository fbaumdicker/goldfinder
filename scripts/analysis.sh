# Installation
mamba create --name coinfinder --file Coinfinder_env.txt
mamba create --name goldfinder --file Goldfinder_env.txt
mamba create --name evo-scope --file EvoScope_env.txt
mamba create --name pandas-basic python pandas matplotlib

git clone git@github.com:fbaumdicker/goldfinder.git
git clone git@github.com:Maxime5G/EvoScope.git
cd EvoScope
make
cd ..

# Run Coinfinder

mamba activate coinfinder
bash coin.sh
mamba deactivate

# Run scripts to create plots that determine appropriate D-value and filter by D- and p-value
mamba activate pandas-basic
python plot_D.py
python filter_d.py  # this has a threshold of D-value == -0.2 hardcoded! Adapt if appropriate
python filter_p.py  # this filters out if p-value > 0.05
mamba deactivate

# Run Goldfinder
mamba activate goldfinder
bash gold.sh  # Uses bonferroni correction - remove resepctive argument to use default of FDR
mamba deactivate

# Run EvoScope
# This did not work for me yet - maybe it is not set up for data this large?
mamba activate EvoScope
bash evo.sh
mamba deactivate
