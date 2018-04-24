#########################
#  Author: Te-Yuan Liu
#########################

#########################
#  Import Library
#########################
rm(list=ls())
library("igraph")

#########################
#  Define Function
#########################
get_dd = function(vec){
    hi = hist(vec, -1:max(vec), plot=FALSE)$count
    hi
}

dd_plot_and_fit = function(graph){
    d = degree(graph, mode="all")
    dd = degree.distribution(graph, mode="all", cumulative=FALSE)
    degree = 1:max(d)
    probability = dd[-1]
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
    print(paste("Alpha = ", round(alpha, 3)))
    print(paste("R square =", round(R.square, 3)))

    plot(probability ~ degree, log="xy", xlab="Degree (log)", ylab="Probability (log)", col=1, main="Degree Distribution")
    curve(power.law.fit, col="red", add=T, n=length(d))
    alpha
}

random_dd_plot_and_fit = function(graph){
    n_v = V(graph)
    n_v_l = length(n_v)
    d = numeric(n_v_l)
    for(i in 1:n_v_l){
        n = n_v[sample(n_v_l, 1)]
        nns = neighbors(graph, n)
        nns_l = length(nns)
        nn = nns[sample(nns_l, 1)]
        d_nn = degree(graph)[nn]
        #print(paste("Step: ", i, ", Node: ", nn, ", Degree: ", d_nn))
        d[i] = d_nn
    }
    dd = get_dd(d)
    degree = 1:max(d)
    probability = dd[-1]
    nonzero.position = which(probability != 0)
    probability = probability[nonzero.position]
    degree = degree[nonzero.position]
    reg = lm(log(probability) ~ log(degree))
    cozf = coef(reg)
    power.law.fit = function(x) exp(cozf[[1]] + cozf[[2]] * log(x))
    alpha = -cozf[[2]]
    R.square = summary(reg)$r.squared
    print(paste("Alpha = ", round(alpha, 3)))
    print(paste("R square =", round(R.square, 3)))
    plot(probability ~ degree, log="xy", xlab="Degree (log)", ylab="Probability (log)", col=1, main="Degree Distribution")
    curve(power.law.fit, col="red", add=T, n=length(d))
    alpha
}

plot_k_i = function(m){
    max_i = 1000
    curve(m*sqrt(max_i/x), from=1, to=1000, col="blue", main="Expected Node Degree versus Node Incoming Timestamp", xlab="Node Incoming Timestamp", ylab="Expected Node Degree", lwd=2)
}

m_routine = function(m_value){
    # Analyze preferential attachment model
    #(a)
    g1 = barabasi.game(1000, m=m_value, directed=F)
    plot(g1, vertex.size=4, vertex.label.cex=0.2)

    #(b)
    fg1 = fastgreedy.community(g1)
    layout1 = layout.fruchterman.reingold(g1)
    plot(fg1, g1, layout=layout1, vertex.size=4, vertex.label.cex=0.2)
    print(modularity(g1, membership(fg1)))

    #(c)
    g2 = barabasi.game(10000, m=m_value, directed=F)
    fg2 = fastgreedy.community(g2)
    print(modularity(g2, membership(fg2)))
    plot(g2, vertex.size=4, vertex.label.cex=0.2)
    layout2 = layout.fruchterman.reingold(g2)
    plot(fg2, g2, layout=layout2, vertex.size=4, vertex.label.cex=0.2)

    #(d)
    dd_plot_and_fit(g1)
    dd_plot_and_fit(g2)

    #(e)
    random_dd_plot_and_fit(g1)
    random_dd_plot_and_fit(g2)

    #(f)
    plot_k_i(m_value)
}

stub_matching = function(){
    g1 = barabasi.game(1000, m=1, directed=F)
    fg1 = fastgreedy.community(g1)
    print(modularity(g1, membership(fg1)))
    layout1 = layout.fruchterman.reingold(g1)
    plot(fg1, g1, layout=layout1, vertex.size=4, vertex.label.cex=0.2)

    d_s = degree(g1)
    n_num = vcount(g1)
    e_num = ecount(g1)
    v = numeric(2*e_num)
    index = 1
    for(i in 1:n_num){
        for(j in 1:d_s[i]){
            v[index] = i
            index = index + 1
        }
    }
    v_p = sample(v)
    v_p_m = matrix(v_p, nrow=e_num, ncol=2)
    g2 = graph_from_edgelist(v_p_m, directed=F)
    g2 = simplify(g2)
    fg2 = fastgreedy.community(g2)
    print(modularity(g2, membership(fg2)))
    layout2 = layout.fruchterman.reingold(g2)
    plot(fg2, g2, layout=layout2, vertex.size=4, vertex.label.cex=0.2)
    dd_plot_and_fit(g1)
    dd_plot_and_fit(g2)
}

power_law_test = function(){
    x = c(3,4,5,6,7)
    x_p = 10^x
    y = numeric(length(x))
    y_index = 0
    for(i in x_p){
        g = barabasi.game(i, m=1, directed=F)
        y_index = y_index + 1
        y[y_index] = dd_plot_and_fit(g)
    }
    plot(y ~ x_p, log="", type="o", xlab="Number of nodes", ylab="Power exponent of power law", col=1, main="Power Law Verification") 
}

#########################
#  Main Function
#########################
main = function(){
    #power_law_test()
    m_routine(1)
    print("*****************************************")
    m_routine(2)
    print("*****************************************")
    m_routine(5)
    print("*****************************************")
    stub_matching()
}
main()


