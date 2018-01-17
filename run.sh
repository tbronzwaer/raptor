make harm CPU=1
echo ''
./RAPTOR model.in dump040 1e19 90 1 1 1
echo ''
python plotter.py
