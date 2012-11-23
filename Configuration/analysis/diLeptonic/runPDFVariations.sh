#!/bin/sh

./load_Analysis -s PDF -f emu_ &
./load_Analysis -s PDF -f ee_ &
./load_Analysis -s PDF -f mumu_ &
wait
echo "Done!"
