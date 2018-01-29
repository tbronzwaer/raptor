make harm CPU=1
echo ' '
./RAPTOR model.in dump040 1e21 60 3 3 1
echo ''
python plotter.py
