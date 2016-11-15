library(data.table)
library(ggplot2)

s = c(1,2,6,0,1)

# change to intervals from 0 to 1, with normalized values
sev <- data.table(
  start = head(seq(0,1,(1/length(s))), -1),
  end = tail(seq(0,1,(1/length(s))), -1),
  value = s / sum(s) # values, normalized
)
setkey(sev, start, end)

#### FIND OVERLAPS

## make into for loop

a <- data.table(
  start = c(0.0, 0.1, 0.8),
  end = c(0.1, 0.8, 1.0),
  value = c(2.0, 5.0, 7.0))
setkey(a, start, end)

b <- data.table(
  start = c(0.0, 0.2, 0.5, 0.6),
  end = c(0.2, 0.5, 0.6, 1.0),
  value = c(3.0, 1.0, 12.0, 1.5))
setkey(b, start, end)

over <- foverlaps(a, b)

over2 <- data.table(
  start = over[, ifelse(start > i.start, start, i.start)],
  end = over[, ifelse(end < i.end, end, i.end)],
  value = over[, value + i.value])

#### PLOTTING
### add a row to over2 for start = 1.0, end = 1.0, value = (as in last row)
p <- ggplot() + geom_step(data = over2, aes(start, value)) + theme(axis.ticks=element_blank(), axis.text.x=element_blank()) + xlab("CDS") + ylab("sum(SHAPE)")

ggsave(file = "/Volumes/USELESS/META/SHAPES/whole_cds/sam24.png")
