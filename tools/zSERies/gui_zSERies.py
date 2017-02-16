'''
GUI tool to convert folder of SER files with EMI to z-stack EMD file.
'''


import sys
import os
import re
import numpy as np
from PyQt4 import QtGui, QtCore

import ncempy.io.ser
import ncempy.io.emd


class Converter(QtGui.QWidget):

    def __init__(self):
        super().__init__()
        
        self.initUI()
        
    def initUI(self):
        #self.statusBar().showMessage('Welcome')
        
        self.resize(600,200)
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())
        
        dir_lbl = QtGui.QLabel('input directory:', self) 
        
        self.dir_txt = QtGui.QLineEdit(self)
        self.dir_txt.setReadOnly(True)
        self.dirButton = QtGui.QPushButton('Open', self)
        
        hbox_dir = QtGui.QHBoxLayout()
        hbox_dir.addWidget(self.dir_txt)
        hbox_dir.addWidget(self.dirButton)       
        
        dir_sep = QtGui.QFrame()
        dir_sep.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        
        emd_lbl = QtGui.QLabel('output EMD file:', self) 
        
        self.emd_txt = QtGui.QLineEdit(self)
        self.emd_txt.setReadOnly(True)
        self.emdButton = QtGui.QPushButton('Save', self)
        
        hbox_emd = QtGui.QHBoxLayout()
        hbox_emd.addWidget(self.emd_txt)
        hbox_emd.addWidget(self.emdButton)
        
        emd_sep = QtGui.QFrame()
        emd_sep.setFrameStyle(QtGui.QFrame.HLine | QtGui.QFrame.Sunken)
        
        
        bin_lbl = QtGui.QLabel('binning factor:', self)
        self.bin_txt = QtGui.QLineEdit(self)
        hbox_bin = QtGui.QHBoxLayout()
        hbox_bin.addWidget(bin_lbl)
        hbox_bin.addWidget(self.bin_txt)
        self.bin_txt.setText('1')
        

        self.msg = QtGui.QLabel('Ready', self)
        self.convButton = QtGui.QPushButton('Convert', self)
        self.exitButton = QtGui.QPushButton('Exit', self)
        
        hbox_buttons = QtGui.QHBoxLayout()
        hbox_buttons.addWidget(self.msg)
        hbox_buttons.addStretch(1)
        hbox_buttons.addWidget(self.convButton)
        hbox_buttons.addWidget(self.exitButton)
  
        
        vbox = QtGui.QVBoxLayout()
        vbox.addWidget(dir_lbl)
        vbox.addLayout(hbox_dir)
        vbox.addWidget(dir_sep)
        vbox.addWidget(emd_lbl)
        vbox.addLayout(hbox_emd)
        vbox.addWidget(emd_sep)
        vbox.addLayout(hbox_bin)
        vbox.addStretch(1)
        vbox.addLayout(hbox_buttons)
        
        self.setLayout(vbox)
        
        self.setWindowTitle('Convert SER to EMD')
        self.show()
        
        
        self.dirButton.clicked.connect(self.clicked_dirButton)
        self.emdButton.clicked.connect(self.clicked_emdButton)
        
        self.convButton.clicked.connect(self.convert)
        self.exitButton.clicked.connect(QtCore.QCoreApplication.instance().quit)
        
        
        
    def keyPressEvent(self, e):
        
        # esc to close
        if e.key() == QtCore.Qt.Key_Escape:
            self.close()
            
    def clicked_dirButton(self):
        self.msg.setText('Ready')
        fname = QtGui.QFileDialog.getExistingDirectory(self, 'Open input folder')
        self.dir_txt.setText(fname)
    
    def clicked_emdButton(self):
        self.msg.setText('Ready')
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save EMD file')
        self.emd_txt.setText(fname)
        
    def convert(self):

        try:
            binfac = int(self.bin_txt.text())
            assert( binfac % 2 == 0 )
        except:
            binfac = 1

        input_dir = self.dir_txt.text()

        rawlist = os.listdir(input_dir)
        
        finelist = []
        for fname in rawlist:
            if fname.endswith('.emi'):
                finelist.append(fname)
                
        assert len(finelist)>0
        
        if os.path.isfile(self.emd_txt.text()):
            os.remove(self.emd_txt.text())
            
        femd = ncempy.io.emd.fileEMD(self.emd_txt.text())
        
        # create EMD group    
        grp = femd.file_hdl['data'].create_group(os.path.basename(input_dir))
        grp.attrs['emd_group_type'] = 1
            
        
        # use first dataset to layout memory
        fser = ncempy.io.ser.fileSER( os.path.join(input_dir, re.sub('\.emi', '_1.ser', finelist[0])), os.path.join(input_dir, finelist[0]) )
        data, first_meta = fser.getDataset(0)
        first_tag = fser.getTag(0)

        dset = grp.create_dataset( 'data', (len(finelist), first_meta['ArrayShape'][1]/binfac, first_meta['ArrayShape'][0]/binfac), dtype=ncempy.io.ser.fileSER.dictDataType[first_meta['DataType']] )
        
        zs = np.zeros(len(finelist))
        
        for i in range(len(finelist)):
            print('importing file {}'.format(finelist[i]))
        
            fser = ncempy.io.ser.fileSER(os.path.join(input_dir, re.sub('\.emi', '_1.ser', finelist[i])), os.path.join(input_dir, finelist[i]) )
        
            assert(fser.head['DataTypeID'] == 0x4122)
            assert(fser.head['ValidNumberElements'] == 1)
          
            data, meta = fser.getDataset(0)
            tag = fser.getTag(0)
            
            assert(meta['DataType'] == first_meta['DataType'])
            
            if not binfac == 1:
                shape = (data.shape[0]/binfac, data.shape[1]/binfac)
                sh = shape[0],data.shape[0]//shape[0],shape[1],data.shape[1]//shape[1]
                data = data.reshape(sh).mean(-1).mean(1)
            
        
            dset[i,:,:] = data[:,:]
            
            zs[i] = fser.emi['Stage Z [um]']
        
        # create dimension datasets
        dims = []
        
        # z
        dims.append( (zs, 'Stage Z', '[um]') )
        
        dim = fser.createDim(first_meta['ArrayShape'][1], first_meta['Calibration'][1]['CalibrationOffset'], first_meta['Calibration'][1]['CalibrationDelta'], first_meta['Calibration'][1]['CalibrationElement'])
        if not binfac ==1:
            shape = (dim.shape[0]/binfac, binfac)
            dim = dim.reshape(shape).mean(-1)
        dims.append( (dim, 'y', '[m]') )
                                                                          
        dim = fser.createDim(first_meta['ArrayShape'][0], first_meta['Calibration'][0]['CalibrationOffset'], first_meta['Calibration'][0]['CalibrationDelta'], first_meta['Calibration'][0]['CalibrationElement'])
        if not binfac ==1:
            shape = (dim.shape[0]/binfac, binfac)
            dim = dim.reshape(shape).mean(-1)
        dims.append( (dim, 'x', '[m]') )

        # write dimensions
        for i in range(len(dims)):
            femd.write_dim('dim{:d}'.format(i+1), dims[i], grp)
                        
        ## write out time as additional dim vector
        #f.write_dim('dim1_time', (time, 'timestamp', '[s]'), grp)
                
        # write comment into Comment group
        femd.put_comment('Combined SER files into a z-stack EMD using the openNCEM tools.')           
        
        self.msg.setText('Done')
            
    
    
        
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    con = Converter()
    sys.exit(app.exec_())
