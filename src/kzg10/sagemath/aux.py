import hashlib

def rho(self, *arg):
    return int(hashlib.sha256((str(arg)).encode()).hexdigest(),16)
