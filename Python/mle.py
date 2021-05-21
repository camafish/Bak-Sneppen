
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

def est_x_k_l(xs):
    lambdas = [(xs[i+1] - xs[i])/(xs[i] - xs[i-1]) for i in range(len(xs) - 1)]

    est_lambda = lambdas[-1]  #sum(lambdas[-2:])/2

    # print(xs)
    # print(lambdas)
    # print(est_lambda)

    kappas = [(xs[i+1] - xs[i])/(est_lambda**i - est_lambda**(i+1)) for i in range(len(xs)-1)]

    est_kappa = kappas[-1]

    # print(kappas)
    # print(est_kappa)

    est_xs = [xs[i] + est_kappa*est_lambda**i for i in range(len(xs))]

    est_x = est_xs[-1]

    print(est_x, est_kappa, est_lambda)

    # plt.scatter(range(len(xs)), est_xs)
    # plt.show()
    return est_x, est_kappa, est_lambda

def model(x,k,l,i):
    return x - k*(l**i)

def ith_term(x,k,l,i, xs, stdevs):
    return (xs[i] - model(x,k,l,i))/(stdevs[i]**2)

def chi2(x,k,l, xs, stdevs):
    return sum([((xs[i] - model(x,k,l,i))**2)/(stdevs[i]**2) for i in range(len(xs))])

def dx(x,k,l, xs, stdevs):
    return -2*sum([ith_term(x,k,l,i, xs, stdevs) for i in range(len(xs))])

def dk(x, k, l,xs,  stdevs):
    return 2*sum([ith_term(x,k,l,i, xs, stdevs)*(l**i) for i in range(len(xs))])

def dl(x, k, l, xs, stdevs):
    return 2*sum([ith_term(x,k,l,i, xs, stdevs)*k*i*(l**(i-1)) for i in range(len(xs))])


# mt_best_fit = [0.6667502352820798 ,0.05765550505778033 , 0.4877514817496227]
# xoshiro_best_fit = [0.6667441393051681, 0.05651643124601967, 0.4913486400654195]

# mt_best_fit = [0.666833530275883, 0.059319217668455285, 0.490726015106669]
# xoshiro_best_fit = [0.6667790635628638, 0.05651631251456804, 0.49134854917019477]

mt_best_fit = [0.6668338497918106, 0.05929183028614936, 0.49081478850426513]
xoshiro_best_fit = [0.6668315508933269, 0.05702934767751706, 0.4956899631270207]


wagners_best_fit = [0.6700251595116361, 0.025248162255475347, 0.6544557120706123]

mt_xs = [0.60653, 0.63885, 0.65309, 0.66033, 0.66347, 0.66512, 0.66598, 0.666371, 0.666681]
mt_stdvs = [0.00916, 0.00265, 0.00144, 0.00074, 0.00025, 0.00013, 0.00005, 0.000039, 0.000032]

xoshiro_xs = [0.60193, 0.63981, 0.6536, 0.66031, 0.66349, 0.66502, 0.66596, 0.666354, 0.666686]
xoshiro_stdev = [0.01197,0.00332,0.00128,0.00033,0.00024,0.00013,0.00006, 0.000024, 0.000024]


# mt_xs = [0.60653, 0.63885, 0.65309, 0.66033, 0.66347, 0.66512, 0.66598, 0.666371]
# mt_stdvs = [0.00916, 0.00265, 0.00144, 0.00074, 0.00025, 0.00013, 0.00005, 0.000039]

# xoshiro_xs = [0.60193, 0.63981, 0.6536, 0.66031, 0.66349, 0.66502, 0.66596, 0.666354]
# xoshiro_stdev = [0.01197,0.00332,0.00128,0.00033,0.00024,0.00013,0.00006, 0.000024]


# wagner's
wagner_xs = [.6402, .6543, .6613, .6618, .6653, .6672, .6681, .6687]
wagner_stdvs = [.0024, .0015, .0006, .0003, .0002, .0001, .0001, .00004]

wagner_omitted_xs = [.6402, .6543, .6653, .6672, .6681, .6687]
wagner_omitted_stdvs = [.0024, .0015, .0002, .0001, .0001, .00004]

# x_init, k_init, l_init =  mt_best_fit[0], mt_best_fit[1], mt_best_fit[2]
# xs = mt_xs
# stdevs = mt_stdvs

#x_init, k_init, l_init =  xoshiro_best_fit[0], xoshiro_best_fit[1], xoshiro_best_fit[2]
#xs = xoshiro_xs
#stdevs = xoshiro_stdev

def test(xs, stdevs, init, miniters):
    x_init, k_init, l_init = init

    maxb = .5
    tau = .5
    #c = .9
    miniters = miniters
    maxiters = miniters*10
    maxdnorm = .001
    maxdiff = 5
    dnorm = maxdnorm
    #print(x_init, k_init, l_init)
    path = []
    x, k, l = x_init, k_init, l_init

    def bdd(ls, bound):
        if abs(max(ls) - min(ls)) < bound:
            return 1
        else:
            return 0

    try:
        j = 0
        delx = dx(x,k,l, xs, stdevs)
        delk = dk(x,k,l, xs, stdevs)
        dell = dl(x,k,l, xs, stdevs)
        norm = np.sqrt(delx**2 + delk**2 + dell**2)

        dnorm = norm
        #while(j < 101 or (not bdd(path[-100:], maxdnorm) and j < maxiters)):
        while(j < miniters + 1 or (norm >= 1 and not bdd(path[-100*10:], maxdiff) and j < maxiters)):
            # if j > 101:
            #     if bdd(path[-100:], maxdnorm):
            #         break
                
            j = j + 1
            
            # given a norm of the gradient
            # d is direction of steepest descent
            d = [-delx/norm, -delk/norm, -dell/norm]

            # print("vals: " +str([x, k, l]))
            # print("dir: "+str(d))
        
            # backtrack line search for step size b
            f1 = chi2(x,k,l,xs,stdevs)
            f2 = lambda t : chi2(x - t*delx, k - t*delk, l - t*dell, xs, stdevs)
            
            b = maxb
            while f2(b) > f1 - (b/2)*(norm**2) :
                b = tau*b
            
            #print("b: " + str(b))

            x = x + b*d[0]
            k = k + b*d[1]
            l = l + b*d[2]

            #compute new gradient here
            oldnorm = norm
            
            delx = dx(x,k,l, xs, stdevs)
            delk = dk(x,k,l, xs, stdevs)
            dell = dl(x,k,l, xs, stdevs)
            norm = np.sqrt(delx**2 + delk**2 + dell**2)

            dnorm = abs(norm - oldnorm)
            if (j % 1000 == 0):
                print(j, maxiters, norm, dnorm)

            #path.append(chi2(x,k,l, xs, stdevs))
            path.append(norm)

    except KeyboardInterrupt:
        pass


    # print(xs)
    # print(x_init, k_init, l_init)
    # print(x, k, l)

    return path, x, k, l

wagner_omitted_init = [
   [0.6698888788465002, 0.009206615856736496, 0.6666667080324012],
   [0.670188500159555, 0.011011398626259699, 0.6661550473912007],
   [0.6701742324053057, 0.011046329994097753, 0.6644242103620056] 
]
# path1, x1, k1, l1 = test(wagner_omitted_xs, wagner_omitted_stdvs, wagner_omitted_init[-1])
# print("wagner: ")
# print("["+str(x1)+", "+str(k1)+", " +str(l1) + "]")

# this is mt
init1 = [
[0.6678479070634551, 0.007611362474209445, 0.7928391734314774],
[0.6684401103548017, 0.010827820456617775, 0.7920691143900779],
[0.6683851857314282, 0.010921104522570117, 0.7875777982051813],
[0.6681677281659393, 0.011797511869998, 0.7662214764014885],
[0.6679318932901442, 0.013303473302129623, 0.7378768799306263],
[0.667813518419899, 0.014400647339750291, 0.7211283311725272],
[0.6678110802809508, 0.014426508895326083, 0.7207623599205217],
[0.6678086322793355, 0.014452313276268036, 0.7203984162182798],
[0.6678062291115638, 0.014478052751358123, 0.720036514670994],
[0.6678038597286193, 0.014503728270923612, 0.7196766407791926],
[0.6678014612532163, 0.014529479171861088, 0.7193169044015594],
[0.6677990987624011, 0.014555166433350702, 0.7189591667409105],
[0.667796784616222, 0.014580782486275435, 0.7186035004249949],
[0.6677741525722264, 0.014835552171436029, 0.7151257608962925],
[0.6675968127296427, 0.017387464945463876, 0.6850601250229323],
[0.667439210030551, 0.02082953163998007, 0.6536902179077224],
[0.6672994236978781, 0.02533147845705583, 0.6218673952866414],
[0.6671796455791131, 0.03076362039346228, 0.5916250473911835],
[0.6670792319588329, 0.03680514014808391, 0.564331885928305],
[0.6669888524756883, 0.04368502583636777, 0.538353712593029],
[0.666980439036306, 0.04440172177855754, 0.5358767496040949],
[0.666931656162162, 0.0488474074941382, 0.5212711477922097],
[0.666879656669248, 0.05414944312743822, 0.5052679612304224],
[0.6668783086660096, 0.0542950453074192, 0.504846641697718],
[0.6668783085496017, 0.054295045308503126, 0.5048466416984654],
[0.6668783086660096, 0.0542950453074192, 0.504846641697718],
[0.6668783085496017, 0.054295045308503126, 0.5048466416984654]
]

# [0.6669391726147968, 0.04813116796978898, 0.523542994892323],
# [0.6669134607929214, 0.05063438739982389, 0.5157242705735483],
# [0.6668967649963208, 0.052338617875073296, 0.5105839458388669],
# [0.666887627386312, 0.053296362637576096, 0.5077544630321752],
# [0.6668827664501168, 0.05381500555356203, 0.5062390624139061],
# [0.6668799104573174, 0.05412184558641407, 0.5053479173143074],
# [0.6668789113390479, 0.0542299683056182, 0.5050348418325641]
# ]

# this is xoshiro
init2 = [
[0.6684539273157956, 0.0070005229741413955, 0.8426397443694927], 
[0.6684457498494315, 0.007009927245716212, 0.842639842607183],
[0.669479940097906, 0.010726790971938604, 0.8413799348328491],
[0.6694238976219828, 0.010763595933617236, 0.8388392306571567],
[0.6693672343594104, 0.010805605225071905, 0.8361664684144594],
[0.6693098309484135, 0.010853309792343793, 0.8333665661637686],
[0.6692538447510485, 0.010904866516688353, 0.8305667276708893],
[0.6691938222188113, 0.010966360690399178, 0.8274606998875835],
[0.6690319132613515, 0.011168974842529105, 0.818482981272035],
[0.6688634338072948, 0.011447473964511816, 0.8081393334289882],
[0.6686466378899313, 0.01193812351521575, 0.7930747166768992],
[0.6685679349615299, 0.012162430286906613, 0.7870553948983999],
[0.6684557814274114, 0.01253576472218502, 0.7778799691182131],
[0.6683103041524299, 0.013133058970403816, 0.7648438541222622],
[0.6681074597468418, 0.01424810879783617, 0.7441923211029695],
[0.6679134298849918, 0.015752521686104196, 0.7212268403605722],
[0.6677290254547296, 0.01777545508513774, 0.6959083432371462],
[0.6675846449870358, 0.019954654626164553, 0.6732797084840011],
[0.667463714581853, 0.022351501757063803, 0.6521672837088854],
[0.6673646862905934, 0.024840427695821437, 0.6332376338623308],
[0.6672845624784925, 0.027306654710069597, 0.6167244356740348],
[0.667216848781513, 0.02978757763593826, 0.6018582938885862],
[0.6671597820183056, 0.03222166605618938, 0.5886340084118524],
[0.667099165928921, 0.03523136930430373, 0.5737902108940188],
[0.6670654074047817, 0.037129202843950625, 0.5651466958705628],
[0.667038126237438, 0.038794631361100215, 0.5579554010674377],
[0.6670107123109361, 0.04060343449599357, 0.5505179018517707],
[0.666994594044332, 0.04173434349280304, 0.5460478925735742],
[0.6669933122543429, 0.04182655413097941, 0.5456891874264501],
[0.6669764466643696, 0.04307330737083778, 0.5409207703147938],
[0.6669423273007526, 0.04579804822818038, 0.5309917546035569],
[0.6669349726043318, 0.04642388913080165, 0.5287992915131047],
[0.6669349726043318, 0.04642388913080165, 0.5287992915131047],
[0.6669349726043318, 0.04642388913080165, 0.5287992915131047]
]

def run_optim(xs, stdvs, init, name, miniters):
    path, x, k, l = test(xs, stdvs, init, miniters)
    print(name)
    print("["+str(x)+", "+str(k)+", " +str(l) + "]")
    plt.plot(range(len(path)), path)
    plt.show()


xoshiro_xs_omit100 = [0.63981, 0.6536, 0.66031, 0.66349, 0.66502, 0.66596, 0.666354, 0.666686]
xoshiro_stdev_omit100 = [0.00332,0.00128,0.00033,0.00024,0.00013,0.00006, 0.000024, 0.000024]

xoshiro_omit_init = [
[0.6669629854126617, 0.02345091392153565, 0.5379956880256407],
[0.6669589383898806, 0.023570667997701633, 0.5368165492586843]
]

# run_optim(mt_xs, mt_stdvs, init1[-1], "MT", 100)
# run_optim(xoshiro_xs_omit100, xoshiro_stdev_omit100, xoshiro_omit_init[-1], "XOSHIRO_OMIT_100", 30*100000)



# to get error bars on x*, take as starting point random normally distributed values within
# mean and stdev as measured.

# wiggled_miniters = 500000

# output = []
# for i in range(10):
#     print("STARTING MT #"+str(i))
#     mt_wiggled = [np.random.normal(mt_xs[i], mt_stdvs[i], 1)[0] for i in range(len(mt_xs))]
#     path, x, k, l = test(mt_wiggled, mt_stdvs, init1[-1], wiggled_miniters)
#     output.append((x, k, l))

# output2 = []
# for i in range(10):
#     print("STARTING XOS #"+str(i))
#     xos_wiggled = [np.random.normal(xoshiro_xs[i], xoshiro_stdev[i], 1)[0] for i in range(len(xoshiro_xs))]
#     path, x, k, l = test(xos_wiggled, xoshiro_stdev, init2[-1], wiggled_miniters)
#     output2.append((x, k, l))

# print()
# print("These are outputs obtained from varying mt xs within their stdevs")
# print(output)

# print("These are outputs obtained from varying xoshiro xs within their stdevs")
# print(output2)



def plot_mle(xs, stdevs, init, name):
    x, k, l = init
    plt.scatter(range(len(xs)), xs)
    interval = np.arange(0, len(xs)-1, .01)
    plt.plot(interval, [model(x, k, l, i) for i in interval], c='red')

    plt.title("Maximum likelihood estimate for threshold model x(k) = x* - CL^k \n x* ~ "+str(np.around(x,5))+", C ~ "+str(np.around(k,5))+", L ~ "+str(np.around(l,5))+"\n"+name)
    plt.xlabel("k (i.e. log_2(N) - log_2(100)")
    plt.ylabel("Threshold")
    plt.errorbar(range(len(xs)), xs, yerr=stdevs, fmt='o')

    plt.show()

wagner_xs = [.6402, .6543, .6613, .6618, .6653, .6672, .6681, .6687]
wagner_stdvs = [.0024, .0015, .0006, .0003, .0002, .0001, .0001, .00004]
wagners_paper_fit = [0.6698, 0.0269, 0.63]

# plot_mle(mt_xs, mt_stdvs, init1[-1], "Mersenne")
plot_mle(xoshiro_xs, xoshiro_stdev, init2[-1], "Xoshiro")
# plot_mle(wagner_xs, wagner_stdvs, wagners_paper_fit, "Wagner")
plot_mle(xoshiro_xs_omit100, xoshiro_stdev_omit100, xoshiro_omit_init[-1], "Xoshiro with n = 100 omitted")


