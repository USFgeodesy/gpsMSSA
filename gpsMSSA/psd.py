from scipy import signal
import matplotlib.pyplot as plt

def psd(x,fs):
    f, Pxx_den = signal.periodogram(x*10, fs)
    plt.semilogy(f, Pxx_den)
    plt.xscale('log')
    plt.ylim(0.0001,1000)
    plt.xlabel('Cycles / year')
    plt.ylabel('PSD [mm**2/cpy]')
    plt.show()

def psd_raw(df,index,sr):
    x = df[df.columns[index]]
    psd(x,sr)
