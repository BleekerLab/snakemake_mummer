# usage: $1 is your reference fasta file, $2 is the query fasta file, $3 is the prefix you want to use
# Align input $1 (reference) to $2 (query)
mummer -mum -n -l 20 -b -c $1 $2 > tmp.mums

# Fix file
sed -i '/^>/! s/^/ /' tmp.mums

# Plot alignment with filename $3
mummerplot --postscript --prefix=$3 tmp.mums

# Generate postscript from plot
gnuplot ${3}.gp

# Convert to pdf
ps2pdf ${3}.ps ${3}.pdf

# Cleanup
rm ${3}.gp ${3}.ps ${3}.fplot ${3}.rplot
