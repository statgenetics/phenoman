#!/usr/bin/env python
#Gao Wang (*) 2012
#This is utility functions extracted from PypeR described by the following document:
#	@article{Xia:McClelland:Wang:2010:JSSOBK:v35c02,
#	  author =	"Xiao-Qin Xia and Michael McClelland and Yipeng Wang",
#	  title =	"PypeR, A Python Package for Using R in Python",
#	  journal =	"Journal of Statistical Software, Code Snippets",
#	  volume =	"35",
#	  number =	"2",
#	  pages =	"1--8",
#	  day =  	"30",
#	  month =	"7",
#	  year = 	"2010",
#	  CODEN =	"JSSOBK",
#	  ISSN = 	"1548-7660",
#	  bibdate =	"2010-03-23",
#	  URL =  	"http://www.jstatsoft.org/v35/c02",
#	  accepted =	"2010-03-23",
#	  acknowledgement = "",
#	  keywords =	"",
#	  submitted =	"2009-10-23",
#	}

import os
import sys
import time
import re
import tempfile
from types import *

try:
	from subprocess import PIPE, Popen
except Exception as e:
	sys.exit(e)

if sys.version < '2.3': # actually python 2.3 is required by tempfile.mkstemp !!!
	set = frozenset = tuple
	basestring = str
elif sys.version < '2.4':
	from sets import Set as set, ImmutableSet as frozenset

if sys.version < '3.0':
	_mystr = _mybytes = lambda s:s
else:
	from functools import reduce
	long, basestring, unicode = int, str, str
	_mybytes = lambda s:bytes(s, 'utf8') #'ascii')
	_mystr = lambda s:str(s, 'utf8')
	StringType = str

try:
	import numpy
	has_numpy = True
except:
	has_numpy = False

_has_subp = False

def NoneStr(obj): return 'NA'

def BoolStr(obj): return obj and 'TRUE' or 'FALSE'

def ReprStr(obj): return repr(obj)

def LongStr(obj):
	rtn = repr(obj)
	if rtn[-1] == 'L': rtn = rtn[:-1]
	return rtn

def ComplexStr(obj):
	return repr(obj).replace('j', 'i')

def SeqStr(obj, head='c(', tail=')'):
	if not obj: return head + tail
	# detect types
	if isinstance(obj, set):
		obj = list(obj)
	# gw: obj_cnone created to remove None objects
	obj_cnone = [x for x in obj if x is not None]
	if len(obj_cnone) == 0:
		return (head + ','.join(list(map(Str4R, obj))) + tail)
	obj0 = obj_cnone[0]
	tp0 = type(obj0)
	simple_types = [str, bool, int, long, float, complex]
	num_types = [int, long, float, complex]
	is_int = tp0 in (int, long) # token for explicit converstion to integer in R since R treat an integer from stdin as double
	if tp0 not in simple_types: head = 'list('
	else:
		tps = isinstance(obj0, basestring) and [StringType] or num_types
		for i in obj_cnone[1:]:
			tp = type(i)
			if tp not in tps:
				head = 'list('
				is_int = False
				break
			elif is_int and tp not in (int, long):
				is_int = False
	# convert
	return (is_int and 'as.integer(' or '') + head + ','.join(list(map(Str4R, obj))) + tail + (is_int and ')' or '')

def DictStr(obj):
	return 'list(' + ','.join(['%s=%s' % (Str4R(a[0]), Str4R(a[1])) for a in list(obj.items())]) + ')'

def OtherStr(obj):
	if has_numpy:
		if isinstance(obj, numpy.ndarray):
			shp = obj.shape
			tpdic = {'i':'as.integer(c(%s))', 'u':'as.integer(c(%s))', 'f':'as.double(c(%s))', 'c':'as.complex(c(%s))', 'b':'c(%s)', 'S':'c(%s)', 'a':'c(%s)', 'U':'c(%s)', 'V':'list(%s)'} # in order: (signed) integer, unsigned integer, float, complex, boolean, string, string, unicode, anything

			def getVec(ary):
				tp = ary.dtype.kind
				rlt = ary.reshape(reduce(lambda a,b=1:a*b, ary.shape))
				rlt = tp == 'b' and [a and 'TRUE' or 'FALSE' for a in rlt] or rlt.tolist()
				if tp != 'V':
					return tpdic.get(tp, 'c(%s)') % repr(rlt)[1:-1]
				# record array
				rlt = list(map(SeqStr, rlt)) # each record will be mapped to vector or list
				return tpdic.get(tp, 'list(%s)') % (', '.join(rlt)) # use str here instead of repr since it has already been converted to str by SeqStr

			if len(shp) == 1: # to vector
				tp = obj.dtype
				if tp.kind != 'V': 
					return getVec(obj)
				# One-dimension record array will be converted to data.frame
				def mapField(f):
					ary = obj[f]
					tp = ary.dtype.kind
					return '"%s"=%s' % (f, tpdic.get(tp, 'list(%s)') % repr(ary.tolist())[1:-1])
				return 'data.frame(%s)' % (', '.join(list(map(mapField, tp.names))))
			elif len(shp) == 2: # two-dimenstion array will be converted to matrix
				return 'matrix(%s, nrow=%d, byrow=TRUE)' % (getVec(obj), shp[0])
			else: # to array
				dim = list(shp[-2:]) # row, col
				dim.extend(shp[-3::-1])
				newaxis = list(range(len(shp)))
				newaxis[-2:] = [len(shp)-1, len(shp)-2]
				return 'array(%s, dim=c(%s))' % (getVec(obj.transpose(newaxis)), repr(dim)[1:-1])
			# record array and char array
	if hasattr(obj, '__iter__'): # for iterators
		if hasattr(obj, '__len__') and len(obj) <= 10000:
			return SeqStr(list(obj))
		else: # waiting for better solution for huge-size containers
			return SeqStr(list(obj))
	return repr(obj)

base_tps = [type(None), bool, int, long, float, complex, str, unicode, list, tuple, set, frozenset, dict] # use type(None) instead of NoneType since the latter cannot be found in the types module in Python 3
base_tps.reverse()
str_func = {type(None):NoneStr, bool:BoolStr, long:LongStr, int:repr, float:repr, complex:ComplexStr, str:repr, unicode:repr, list:SeqStr, tuple:SeqStr, set:SeqStr, frozenset:SeqStr, dict:DictStr}

def Str4R(obj):
	'''
	convert a Python basic object into an R object in the form of string.
	'''
	#return str_func.get(type(obj), OtherStr)(obj)
	if type(obj) in str_func:
		return str_func[type(obj)](obj)
	for tp in base_tps:
		if isinstance(obj, tp):
			return str_func[tp](obj)
	return OtherStr(obj)

def runR(maer):
	try:
		warnings = None
		if not isinstance(maer, basestring):
			raise ValueError("Input script must be a string")
		tc = Popen(["R", '--subordinate', '--vanilla'], stdin = PIPE, stdout = PIPE, stderr = PIPE)
		out, error = tc.communicate(maer.encode(sys.getdefaultencoding()))
		if (tc.returncode):
			raise ValueError("**** R ERROR ****\n{}*****************\n".format(error.decode(sys.getdefaultencoding())))
		else:
			if error:
				warnings = "**** R WARNING ****\n{}*******************\n".format(error.decode(sys.getdefaultencoding()))
	except Exception as e:
		sys.stdout.write(str(e))
		raise ValueError('See **** R ERROR **** above!')
	return out.decode(sys.getdefaultencoding()), warnings
