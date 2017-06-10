revenue <- c(30,20, 23, 17)
product <- factor(c("bread", "cake", "bread", "cake"))
shop <- gl(2,2, labels=c("shop_1", "shop_2"))

(shop_revenue <- ave(revenue, shop, FUN=sum))
(shop_revenue2 <- ave(revenue, product, FUN=sum))


ave

interaction(product)
lapply(split(revenue,product), sum)
