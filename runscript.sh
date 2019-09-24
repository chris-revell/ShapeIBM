
for i in $(seq 0.0 0.1 1.0); do
  for j in $(seq 0.5 0.1 1.5); do
    echo $i 
    echo $j
    ./shapeIBM 256 1.0 1.0 0.75 0.2 $i 0.5 50.0 10.0 20 $j 1 0.1 1500 100 1 0 0 >> /dev/null
  done
done

