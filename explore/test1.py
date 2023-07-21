from scipy import stats
import matplotlib.pyplot as plt
def Campbell(p, np,wp, bo, cf, cw, sw0, boi, date,pi):
    Eo = (bo - boi)
    Efw = boi * (((cw * sw0) + cf) / (1 - sw0)) * (pi - p)
    F = (np * bo) + (wp)
    y = F/(Eo+Efw)
    x = date
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y)
    # Nombres de los ejes
    plt.title('Gráfico Campbell')
    plt.xlabel("Date (Year) ")
    plt.ylabel('F/Eo+Efw ')
    plt.show(
    )

def G_method(pr, np,wp, bo, cf, cw, sw0, boi,pi,we):
    Eo = (bo - boi)
    Efw = boi * (((cw * sw0) + cf) / (1 - sw0)) * (pi - pr)
    F = (np * bo) + (wp)
    y = F / (Eo + Efw)
    x = we/ (Eo + Efw)
    slope, intercept, r, p, se = stats.linregress(x, y)
    N1 = slope
    yt = intercept + (x * N1)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(x, y)
    plt.plot(x, y)
    ax.plot(x, yt, c='g', label='Regresion lineal')
    # Nombres de los ejes
    text = " N [MMStb]: %.3f " % (intercept / 1000000)
    plt.title('Gráfico Havlena y Odeh' )
    plt.xlabel("WeBw/Eo+Efw")
    plt.ylabel('F/Eo+Efw ')
    plt.text(1e+7, 2e+8, text)
    plt.legend()
    plt.show()