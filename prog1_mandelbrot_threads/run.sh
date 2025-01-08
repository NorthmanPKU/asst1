for i in {2..16}; do
  speedup=$(./mandelbrot -t $i | grep -oP '\(\K[0-9.]+(?=x speedup)')
  echo "$i $speedup" >> result
done