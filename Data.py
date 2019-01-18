class Data():
    def __init__(self,filename,n_spins,q):
        self.data = np.genfromtxt(filename)
        self.n_samples = self.data.shape[0]
        self.n_spins = n_spins
        self.n_states = n_states
        assert n_spins*n_states == self.data[0].size
        self.data = self.data.reshape(self.n_samples,L,q)
        self.split(0.8,0.1,1)
    def split(p_train,p_test,seed):
        r = np.random.rand(size=self.n_samples,seed=seed)
        self.ix_train = np.where(r < p_train)
        self.ix_val = np.where(r > p_train and r < p_train + p_test)
        self.ix_test = np.where(r > p_train + p_test)
