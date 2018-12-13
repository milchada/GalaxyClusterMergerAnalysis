from astropy.io import fits

class Observation:
	def __init__(self, name):
		self.name = name

	def add_data(datafile):
		self.data = fits.getdata(datafile)

	def add_error(errorfile):
		self.error = fits.getdata(errorfile)