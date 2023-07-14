# Load the required packages
require(ggplot2)
require(ggseqlogo)

# Some sample data
data(ggseqlogo_sample)

seqs_dna
ggplot() + geom_logo( seqs_dna$MA0001.1 ) + theme_void()
ggsave("logo-tflogo.png", bg="transparent", width=2, height=1.5)
