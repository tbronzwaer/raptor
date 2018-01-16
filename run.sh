make harm CPU=1
echo ''
./RAPTOR model.in dump040 1e19 60 1 1 1
echo ''
python plotter.py
