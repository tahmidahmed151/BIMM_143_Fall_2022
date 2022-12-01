Class 05: Data Visualization with GGPLOT
================
Tahmid Ahmed

``` r
#install.packages("ggplot2")
library(ggplot2)
ggplot(cars) + aes(x=speed, y=dist)
```

![](class05_files/figure-gfm/unnamed-chunk-1-1.png)

``` r
ggplot(mpg) + aes(x=displ, y=hwy) +geom_point()
```

![](class05_files/figure-gfm/unnamed-chunk-2-1.png)

``` r
ggplot(mpg) + aes(x=displ, y=hwy) +geom_point() + geom_smooth(method= lm, se= FALSE)
```

    `geom_smooth()` using formula = 'y ~ x'

![](class05_files/figure-gfm/unnamed-chunk-3-1.png)

``` r
ggplot(cars) + aes(x=speed, y=dist) + geom_point() + geom_smooth(method=lm, se= FALSE)
```

    `geom_smooth()` using formula = 'y ~ x'

![](class05_files/figure-gfm/unnamed-chunk-4-1.png)

``` r
url <- "https://bioboot.github.io/bimm143_S20/class-material/up_down_expression.txt"
genes <- read.delim(url)
head(genes)
```

            Gene Condition1 Condition2      State
    1      A4GNT -3.6808610 -3.4401355 unchanging
    2       AAAS  4.5479580  4.3864126 unchanging
    3      AASDH  3.7190695  3.4787276 unchanging
    4       AATF  5.0784720  5.0151916 unchanging
    5       AATK  0.4711421  0.5598642 unchanging
    6 AB015752.4 -3.6808610 -3.5921390 unchanging

``` r
ggplot(genes) + aes(x=Condition1, y=Condition2) + geom_point()
```

![](class05_files/figure-gfm/unnamed-chunk-6-1.png)

``` r
ggplot(genes) + aes(x=Condition1, y=Condition2, col=State) + geom_point() 
```

![](class05_files/figure-gfm/unnamed-chunk-7-1.png)

``` r
 ggplot(genes) + aes(x=Condition1, y=Condition2, col=State) + geom_point() + scale_colour_manual( values=c("blue","gray","red") )
```

![](class05_files/figure-gfm/unnamed-chunk-8-1.png)

``` r
ggplot(genes) + aes(x=Condition1, y=Condition2, col=State) + geom_point() + scale_colour_manual(values=c("blue","gray","red")) +
    labs(title="Gene Expresion Changes Upon Drug Treatment",
         x="Control (no drug) ",
         y="Drug Treatment")
```

![](class05_files/figure-gfm/unnamed-chunk-9-1.png)
