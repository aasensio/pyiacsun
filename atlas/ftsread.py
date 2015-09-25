#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Author: cdiazbas@iac.es
# Date: 16.04.2015
# Code: Translation of IDL ftsread


def ftsread(ini,endi):

    import os
    #import os.path
    from struct import unpack
    import numpy

    #Mensaje importante
    print('Wavelength range (3290 - 12508 A)')
    
    # Directorio actual
    currdir = os.getcwd()

    #Directorio del atlas
    sdir = '/usr/pkg/rsi/idl_local/data'

    #Directorio de trabajo
    tdir = currdir

    #Copiamos los datos necesarios
    skip= (ini-3290)
    count = abs(ini-endi)
    os.system('dd if='+sdir+'/fts_cent.dat of='+tdir+'/tmp bs=1000 skip='+str(skip)+' count='+str(count))

    #Accedemos a los datos
    FileTmp = open('tmp', 'rb')
    FileRead =FileTmp.read()
    VarDecode=unpack('!'+str(len(FileRead)/2)+'h',FileRead)
    VarFinal = numpy.array(VarDecode)

    #Eliminamos el fichero al final
    os.remove(tdir+'/tmp')

    return VarFinal



if __name__ == '__main__':
    
    import pylab
    Atlas = ftsread(ini= 10825,endi = 10843)
        
    pylab.plot(Atlas/1e4)
    pylab.ylim(0,1)
    pylab.show()
    