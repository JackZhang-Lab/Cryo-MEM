# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MembraneAnalyzer_Bezierfit.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MembraneAnalyzer_Bezierfit(object):
    def setupUi(self, MembraneAnalyzer_Bezierfit):
        MembraneAnalyzer_Bezierfit.setObjectName("MembraneAnalyzer_Bezierfit")
        MembraneAnalyzer_Bezierfit.resize(441, 611)
        self.gridLayoutWidget = QtWidgets.QWidget(MembraneAnalyzer_Bezierfit)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(40, 10, 361, 103))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout.setContentsMargins(0, 0, 0, 0)
        self.gridLayout.setObjectName("gridLayout")
        self.particle_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.particle_lineEdit.setObjectName("particle_lineEdit")
        self.gridLayout.addWidget(self.particle_lineEdit, 0, 1, 1, 1)
        self.output_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.output_lineEdit.setObjectName("output_lineEdit")
        self.gridLayout.addWidget(self.output_lineEdit, 2, 1, 1, 1)
        self.template_label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.template_label.setObjectName("template_label")
        self.gridLayout.addWidget(self.template_label, 1, 0, 1, 1)
        self.template_browse_pushButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.template_browse_pushButton.setObjectName("template_browse_pushButton")
        self.gridLayout.addWidget(self.template_browse_pushButton, 1, 2, 1, 1)
        self.output_label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.output_label.setObjectName("output_label")
        self.gridLayout.addWidget(self.output_label, 2, 0, 1, 1)
        self.particle_browse_pushButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.particle_browse_pushButton.setObjectName("particle_browse_pushButton")
        self.gridLayout.addWidget(self.particle_browse_pushButton, 0, 2, 1, 1)
        self.template_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.template_lineEdit.setObjectName("template_lineEdit")
        self.gridLayout.addWidget(self.template_lineEdit, 1, 1, 1, 1)
        self.particle = QtWidgets.QLabel(self.gridLayoutWidget)
        self.particle.setObjectName("particle")
        self.gridLayout.addWidget(self.particle, 0, 0, 1, 1)
        self.output_pushButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.output_pushButton.setObjectName("output_pushButton")
        self.gridLayout.addWidget(self.output_pushButton, 2, 2, 1, 1)
        self.line = QtWidgets.QFrame(MembraneAnalyzer_Bezierfit)
        self.line.setGeometry(QtCore.QRect(40, 180, 361, 16))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.gridLayoutWidget_2 = QtWidgets.QWidget(MembraneAnalyzer_Bezierfit)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(40, 190, 361, 85))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.coarsefit_points_number_label = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.coarsefit_points_number_label.setObjectName("coarsefit_points_number_label")
        self.gridLayout_2.addWidget(self.coarsefit_points_number_label, 0, 0, 1, 1)
        self.coarsefit_iter_label = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.coarsefit_iter_label.setObjectName("coarsefit_iter_label")
        self.gridLayout_2.addWidget(self.coarsefit_iter_label, 1, 0, 1, 1)
        self.coarsefit_points_number_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.coarsefit_points_number_lineEdit.setObjectName("coarsefit_points_number_lineEdit")
        self.gridLayout_2.addWidget(self.coarsefit_points_number_lineEdit, 0, 1, 1, 1)
        self.coarsefit_iter_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.coarsefit_iter_lineEdit.setObjectName("coarsefit_iter_lineEdit")
        self.gridLayout_2.addWidget(self.coarsefit_iter_lineEdit, 1, 1, 1, 1)
        self.coarsefit_cpus_label = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.coarsefit_cpus_label.setObjectName("coarsefit_cpus_label")
        self.gridLayout_2.addWidget(self.coarsefit_cpus_label, 2, 0, 1, 1)
        self.coarsefit_cpus_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_2)
        self.coarsefit_cpus_lineEdit.setObjectName("coarsefit_cpus_lineEdit")
        self.gridLayout_2.addWidget(self.coarsefit_cpus_lineEdit, 2, 1, 1, 1)
        self.gridLayoutWidget_3 = QtWidgets.QWidget(MembraneAnalyzer_Bezierfit)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(40, 120, 361, 61))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.gridLayout_3.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.physical_membrane_distance_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.physical_membrane_distance_lineEdit.setObjectName("physical_membrane_distance_lineEdit")
        self.gridLayout_3.addWidget(self.physical_membrane_distance_lineEdit, 1, 1, 1, 1)
        self.bezier_curve_degree_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_3)
        self.bezier_curve_degree_lineEdit.setObjectName("bezier_curve_degree_lineEdit")
        self.gridLayout_3.addWidget(self.bezier_curve_degree_lineEdit, 0, 1, 1, 1)
        self.bezier_curve_degree_label = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.bezier_curve_degree_label.setObjectName("bezier_curve_degree_label")
        self.gridLayout_3.addWidget(self.bezier_curve_degree_label, 0, 0, 1, 1)
        self.physical_membrane_distance_label = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.physical_membrane_distance_label.setObjectName("physical_membrane_distance_label")
        self.gridLayout_3.addWidget(self.physical_membrane_distance_label, 1, 0, 1, 1)
        self.line_2 = QtWidgets.QFrame(MembraneAnalyzer_Bezierfit)
        self.line_2.setGeometry(QtCore.QRect(40, 270, 361, 16))
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.gridLayoutWidget_4 = QtWidgets.QWidget(MembraneAnalyzer_Bezierfit)
        self.gridLayoutWidget_4.setGeometry(QtCore.QRect(40, 280, 361, 121))
        self.gridLayoutWidget_4.setObjectName("gridLayoutWidget_4")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.gridLayoutWidget_4)
        self.gridLayout_4.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.curve_penalty_thr_label = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.curve_penalty_thr_label.setObjectName("curve_penalty_thr_label")
        self.gridLayout_4.addWidget(self.curve_penalty_thr_label, 0, 0, 1, 1)
        self.dithering_range_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_4)
        self.dithering_range_lineEdit.setObjectName("dithering_range_lineEdit")
        self.gridLayout_4.addWidget(self.dithering_range_lineEdit, 1, 1, 1, 1)
        self.refine_iter_label = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.refine_iter_label.setObjectName("refine_iter_label")
        self.gridLayout_4.addWidget(self.refine_iter_label, 2, 0, 1, 1)
        self.curve_penalty_thr_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_4)
        self.curve_penalty_thr_lineEdit.setObjectName("curve_penalty_thr_lineEdit")
        self.gridLayout_4.addWidget(self.curve_penalty_thr_lineEdit, 0, 1, 1, 1)
        self.dithering_range_label = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.dithering_range_label.setObjectName("dithering_range_label")
        self.gridLayout_4.addWidget(self.dithering_range_label, 1, 0, 1, 1)
        self.refine_iter_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_4)
        self.refine_iter_lineEdit.setObjectName("refine_iter_lineEdit")
        self.gridLayout_4.addWidget(self.refine_iter_lineEdit, 2, 1, 1, 1)
        self.refine_cpus_label = QtWidgets.QLabel(self.gridLayoutWidget_4)
        self.refine_cpus_label.setObjectName("refine_cpus_label")
        self.gridLayout_4.addWidget(self.refine_cpus_label, 3, 0, 1, 1)
        self.refine_cpus_lineEdit = QtWidgets.QLineEdit(self.gridLayoutWidget_4)
        self.refine_cpus_lineEdit.setObjectName("refine_cpus_lineEdit")
        self.gridLayout_4.addWidget(self.refine_cpus_lineEdit, 3, 1, 1, 1)
        self.line_3 = QtWidgets.QFrame(MembraneAnalyzer_Bezierfit)
        self.line_3.setGeometry(QtCore.QRect(40, 430, 361, 16))
        self.line_3.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_3.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_3.setObjectName("line_3")
        self.LOG_label = QtWidgets.QLabel(MembraneAnalyzer_Bezierfit)
        self.LOG_label.setGeometry(QtCore.QRect(210, 440, 31, 16))
        self.LOG_label.setObjectName("LOG_label")
        self.LOG_textBrowser = QtWidgets.QTextBrowser(MembraneAnalyzer_Bezierfit)
        self.LOG_textBrowser.setGeometry(QtCore.QRect(10, 460, 421, 141))
        self.LOG_textBrowser.setObjectName("LOG_textBrowser")
        self.horizontalLayoutWidget = QtWidgets.QWidget(MembraneAnalyzer_Bezierfit)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(40, 400, 361, 32))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.kill_pushButton = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.kill_pushButton.setObjectName("kill_pushButton")
        self.horizontalLayout.addWidget(self.kill_pushButton)
        self.launch_pushButton = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.launch_pushButton.setObjectName("launch_pushButton")
        self.horizontalLayout.addWidget(self.launch_pushButton)

        self.retranslateUi(MembraneAnalyzer_Bezierfit)
        QtCore.QMetaObject.connectSlotsByName(MembraneAnalyzer_Bezierfit)

    def retranslateUi(self, MembraneAnalyzer_Bezierfit):
        _translate = QtCore.QCoreApplication.translate
        MembraneAnalyzer_Bezierfit.setWindowTitle(_translate("MembraneAnalyzer_Bezierfit", "Membrane Analyzer - Bezierfit - MemXTerminator"))
        self.template_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Template .cs file"))
        self.template_browse_pushButton.setText(_translate("MembraneAnalyzer_Bezierfit", "Browse..."))
        self.output_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Output JSON file"))
        self.particle_browse_pushButton.setText(_translate("MembraneAnalyzer_Bezierfit", "Browse..."))
        self.particle.setText(_translate("MembraneAnalyzer_Bezierfit", "Particle .cs file"))
        self.output_pushButton.setText(_translate("MembraneAnalyzer_Bezierfit", "Browse..."))
        self.coarsefit_points_number_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Coarsefit points number"))
        self.coarsefit_iter_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Coarsefit iterations"))
        self.coarsefit_points_number_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "600"))
        self.coarsefit_iter_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "300"))
        self.coarsefit_cpus_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Coarsefit cpus"))
        self.coarsefit_cpus_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "20"))
        self.physical_membrane_distance_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "35"))
        self.bezier_curve_degree_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "3"))
        self.bezier_curve_degree_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Beizer curve degree"))
        self.physical_membrane_distance_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Physical membrane distance(Å)"))
        self.curve_penalty_thr_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Curve penalty threshold"))
        self.dithering_range_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "50"))
        self.refine_iter_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Refine iterations"))
        self.curve_penalty_thr_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "0.05"))
        self.dithering_range_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Dithering range"))
        self.refine_iter_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "700"))
        self.refine_cpus_label.setText(_translate("MembraneAnalyzer_Bezierfit", "Refine cpus"))
        self.refine_cpus_lineEdit.setText(_translate("MembraneAnalyzer_Bezierfit", "12"))
        self.LOG_label.setText(_translate("MembraneAnalyzer_Bezierfit", "LOG"))
        self.kill_pushButton.setText(_translate("MembraneAnalyzer_Bezierfit", "Kill"))
        self.launch_pushButton.setText(_translate("MembraneAnalyzer_Bezierfit", "Launch"))
