from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QMessageBox, QApplication, QDialog
from memxterminator.GUI.radon_gui import Ui_Form
from memxterminator.GUI.radonfit_membrane_analyzer_gui import Ui_MembraneAnalyzer
from memxterminator.GUI.radonfit_membrane_subtraction_gui import Ui_MembraneSubtraction
from memxterminator.GUI.MicrographMembraneSubtraction_Radonfit import Ui_MicrographMembraneSubtraction_Radonfit
import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
import os
import json
import mrcfile
from memxterminator.radonfit.lib._utils import *
from memxterminator.radonfit.lib.radonanalyser import *
import subprocess


class RadonApp(QtWidgets.QDialog, Ui_Form):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        
        self.mrcPathLineEdit.setToolTip("Path to the 2D average MRC file you want to analyze.")
        self.sectionSpinBox.setToolTip("Select the section number you want to analyze from the MRC file.")
        self.cropRateLineEdit.setToolTip("Enter the crop rate value (between 0 and 1). This determines how much of the image is considered for analysis.")
        self.thrLineEdit.setToolTip("Enter the threshold value (between 0 and 1). This determines the sensitivity of the analysis.")
        self.jsonPathLineEdit.setToolTip("Path where the radonanalysis_info.json file will be saved. Default is the directory of the MRC file.")
        
        self.browseButton.clicked.connect(self.browse_file)
        self.previewButton.clicked.connect(self.preview_section)
        self.analyzeButton.clicked.connect(self.analyze_section)
        self.saveButton.clicked.connect(self.save_results)
        self.jsonBrowseButton.clicked.connect(self.browse_json_save_path)
        self.mrc_file = None
        self.image = None
        self.analyzer = None
        self.params = {}

    def browse_file(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open MRC File", "", "MRC Files (*.mrc)")
        if filepath:
            self.mrcPathLineEdit.setText(filepath)
            self.mrc_file = filepath
            total_sections = mrcfile.open(self.mrc_file).data.shape[0]
            self.sectionSpinBox.setMaximum(total_sections - 1)
            for i in range(total_sections):
                self.params[str(i)] = {
                    'crop_rate': 0.6,
                    'thr': 0.7, 
                    'theta_start': 0,
                    'theta_end': 180
                }

    def preview_section(self):
        if self.mrc_file:
            section = self.sectionSpinBox.value()
            self.image = readmrc(self.mrc_file, section=section, mode='gpu')
            plt.imshow(self.image.get(), cmap='gray')
            plt.show()

    def analyze_section(self):
        try:
            if self.image is None and self.mrc_file:
                section = self.sectionSpinBox.value()
                self.image = readmrc(self.mrc_file, section=section, mode='gpu')
            
            if self.image is not None:
                crop_rate = float(self.cropRateLineEdit.text())
                thr = float(self.thrLineEdit.text())
                theta_start = int(self.theta_start_lineEdit.text())
                theta_end = int(self.theta_end_lineEdit.text())
                self.analyzer = RadonAnalyzer('None', 0, self.image, crop_rate, thr, theta_start=theta_start, theta_end=theta_end)
                self.analyzer.visualize_analyze()
                section = self.sectionSpinBox.value()
                self.params[str(section)] = {
                    'crop_rate': crop_rate,
                    'thr': thr, 
                    'theta_start': theta_start,
                    'theta_end': theta_end
                }
        except Exception as e:
            print(f"Error: {e}")
            QMessageBox.warning(self, "Error", "Radon analysis failed. Please try again with different parameters.")
    def browse_json_save_path(self):
        default_dir = QtCore.QFileInfo(self.mrcPathLineEdit.text()).absolutePath()
        save_path, _ = QFileDialog.getSaveFileName(self, "Select JSON Save Path", default_dir, "JSON Files (*.json);;All Files (*)")
        if save_path:
            self.jsonPathLineEdit.setText(save_path)
    def save_results(self):
        save_path = self.jsonPathLineEdit.text() or 'radonanalysis_info.json'
        if os.path.exists(save_path):
            with open(save_path, 'r') as file:
                data = json.load(file)
        else:
            data = {}
        data.update(self.params)
        with open(save_path, 'w') as file:
            json.dump(data, file, indent=4)
        QMessageBox.information(self, "Info", "Results saved successfully!")

class MembraneAnalyzerApp(QtWidgets.QDialog, Ui_MembraneAnalyzer):
    PID_FILE = 'process.pid'
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.templates_starfile_browse_button.clicked.connect(self.browse_templates_starfile)
        self.particles_starfile_browse_button.clicked.connect(self.browse_particles_starfile)
        self.output_starfile_browse_button.clicked.connect(self.output_starfile)
        self.info_JSON_file_browse_button.clicked.connect(self.info_json_file)
        self.kappa_templates_checkBox.stateChanged.connect(self.generate_kappa_templates)


        self.templates_starfile_name = None
        self.particles_starfile_name = None
        self.output_starfile_name = None
        self.info_json_name = None
        self.kappa_template = False
        self.process = None
        self.kappa_number = self.kappa_num_textedit.text()
        self.kappa_start_value = self.kappa_start_textedit.text()
        self.kappa_end_value = self.kappa_end_textedit.text()

        self.initialsigma1 = self.initialsigma1_textedit.text()
        self.initialsigma2 = self.initialsigma2_textedit.text()
        self.template_size = self.templatesize.text()
        self.sigmarange = self.select_range_spinbox.value()
        self.sigma_step = self.sigma_step_for_centerfit_textedit.text()
        self.curve_kappa_start = self.curve_kappa_start_textedit.text()
        self.curve_kappa_end = self.curve_kappa_end_textedit.text()
        self.curve_kappa_step = self.curve_kappa_step_textedit.text()
        self.edge_sigma_for_mask = self.edge_sigma_for_mask_textedit.text()
        self.extra_mem_dist = self.extra_mem_dist_textedit.text()
        self.mem_edge_sigma = self.edge_sigma_for_mem_average_textedit.text()

        self.check_running_process()
        self.Launch_button.clicked.connect(self.start_process)
        self.Kill_button.clicked.connect(self.kill_process)
        
        with open('run.out', 'a') as f:
            f.write(f"Membrane Analysis ready.\n")
        self.last_read_position = 0
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_log)
        self.timer.start(1000)
        

    def browse_templates_starfile(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open STAR File", "", "STAR Files (*.star)")
        if filepath:
            self.templates_starfile_path_textedit.setText(filepath)
            self.templates_starfile_name = filepath

    def browse_particles_starfile(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open STAR File", "", "STAR Files (*.star)")
        if filepath:
            self.particles_starfile_path_textedit.setText(filepath)
            self.particles_starfile_name = filepath
    def output_starfile(self):
        filepath, _ = QFileDialog.getSaveFileName(self, "Save STAR File", "", "STAR Files (*.star)")
        if filepath:
            self.output_starfile_path_textedit.setText(filepath)
            self.output_starfile_name = filepath
    def info_json_file(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open JSON File", "", "JSON Files (*.json)")
        if filepath:
            self.info_json_file_textedit.setText(filepath)
            self.info_json_name = filepath
    def generate_kappa_templates(self, state):
        if state == QtCore.Qt.Checked:
            self.kappa_template = self.whichtemplate_spinBox.value()
        elif state == QtCore.Qt.Unchecked:
            self.kappa_template = False
    def start_process(self):
        self.kappa_number = self.kappa_num_textedit.text()
        self.kappa_start_value = self.kappa_start_textedit.text()
        self.kappa_end_value = self.kappa_end_textedit.text()

        self.initialsigma1 = self.initialsigma1_textedit.text()
        self.initialsigma2 = self.initialsigma2_textedit.text()
        self.template_size = self.templatesize.text()
        self.sigmarange = self.select_range_spinbox.value()
        self.sigma_step = self.sigma_step_for_centerfit_textedit.text()
        self.curve_kappa_start = self.curve_kappa_start_textedit.text()
        self.curve_kappa_end = self.curve_kappa_end_textedit.text()
        self.curve_kappa_step = self.curve_kappa_step_textedit.text()
        self.edge_sigma_for_mask = self.edge_sigma_for_mask_textedit.text()
        self.extra_mem_dist = self.extra_mem_dist_textedit.text()
        self.mem_edge_sigma = self.edge_sigma_for_mem_average_textedit.text()
        params = ['--templates_starfile_name', f'{self.templates_starfile_name}',
                  '--output_filename', f'{self.output_starfile_name}',
                  '--particle_starfile_name', f'{self.particles_starfile_name}',
                  '--kappa_template', f'{int(self.kappa_template)}',
                  '--kappanum', f'{self.kappa_number}',
                  '--kappastart', f'{self.kappa_start_value}',
                  '--kappaend', f'{self.kappa_end_value}',
                  '--info_json', f'{self.info_json_name}',
                  '--sigma1', f'{self.initialsigma1}',
                  '--sigma2', f'{self.initialsigma2}',
                  '--template_size', f'{self.template_size}',
                  '--sigma_range', f'{self.sigmarange}',
                  '--sigma_step', f'{self.sigma_step}',
                  '--curve_kappa_start', f'{self.curve_kappa_start}',
                  '--curve_kappa_end', f'{self.curve_kappa_end}',
                  '--curve_kappa_step', f'{self.curve_kappa_step}',
                  '--edge_sigma', f'{self.edge_sigma_for_mask}',
                  '--extra_mem_dist', f'{self.extra_mem_dist}',
                  '--mem_edge_sigma', f'{self.mem_edge_sigma}']
        # self.process = subprocess.Popen(['python', 'membrane_analysis-main.py'] + params)

        with open('run.out', 'w') as f:
            # self.process = subprocess.Popen(['python','-u', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'membrane_analysis-main.py')] + params, stdout=f, stderr=subprocess.STDOUT)
            self.process = subprocess.Popen(['python','-u', '-m', 'memxterminator.radonfit.bin.membrane_analysis-main'] + params, stdout=f, stderr=subprocess.STDOUT)
        print("Membrane Analysis started with PID:", self.process.pid)
        with open(self.PID_FILE, 'w') as f:
            f.write(str(self.process.pid))
    def kill_process(self):
        if self.process:
            self.process.terminate()
            print(f"Process PID {self.process.pid} terminated")
            self.process = None
            if os.path.exists(self.PID_FILE):
                os.remove(self.PID_FILE)
            self.timer.stop()
    def update_log(self):
        # 读取日志文件内容
        try:
            with open('run.out', 'r') as f:
                f.seek(self.last_read_position)  # 跳转到上次读取的位置
                new_content = f.read()  # 读取新内容
                self.last_read_position = f.tell()  # 更新读取的位置
            if new_content:
                self.textBrowser_log.append(new_content)
        except FileNotFoundError:
            self.textBrowser_log.append("Error: 'run.out' file not found.")
    def check_running_process(self):
        if os.path.exists(self.PID_FILE):
            with open(self.PID_FILE, 'r') as f:
                pid = int(f.read().strip())
            try:
                os.kill(pid, 0)  # Check if process is running
                self.process = pid
                print(f"Process with PID {pid} is still running!")
                # Here, you can also update the GUI to show that the process is running
            except OSError:
                print(f"Process with PID {pid} is not running.")
                os.remove(self.PID_FILE)

class MembraneSubtractionApp(QtWidgets.QDialog, Ui_MembraneSubtraction):
    PID_FILE = 'process.pid'
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.particles_selected_starfile_button.clicked.connect(self.browse_particles_starfile)
        self.membrane_analysis_results_button.clicked.connect(self.browse_mem_analysis_starfile)


        self.mem_analysis_starfile_name = None
        self.particles_starfile_name = None
        self.bias = self.Bias_lineEdit.text()
        self.extra_mem_dist = self.Extra_mem_dist_lineEdit.text()
        self.scaling_factor_start = self.Scaling_factor_start_lineEdit.text()
        self.scaling_factor_end = self.Scaling_factor_end_lineEdit.text()
        self.scaling_factor_step = self.Step_lineEdit.text()
        self.cpu = self.CPU_lineEdit.text()
        self.batch_size = self.Batch_size_lineEdit.text()
        self.process = None

        self.check_running_process()
        self.launch_button.clicked.connect(self.start_process)
        self.kill_button.clicked.connect(self.kill_process)
        
        with open('run.out', 'a') as f:
            f.write(f"Membrane Subtraction ready.\n")
        self.last_read_position = 0
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_log)
        self.timer.start(1000)  # 每秒更新一次
        

    def browse_mem_analysis_starfile(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open STAR File", "", "STAR Files (*.star)")
        if filepath:
            self.membrane_analysis_file_lineEdit.setText(filepath)
            self.mem_analysis_starfile_name = filepath

    def browse_particles_starfile(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open STAR File", "", "STAR Files (*.star)")
        if filepath:
            self.particles_selected_starfile_lineEdit.setText(filepath)
            self.particles_starfile_name = filepath

    def start_process(self):
        self.mem_analysis_starfile_name = self.membrane_analysis_file_lineEdit.text()
        self.particles_starfile_name = self.particles_selected_starfile_lineEdit.text()
        params = ['--particles_selected_filename', f'{self.particles_starfile_name}',
                  '--membrane_analysis_filename', f'{self.mem_analysis_starfile_name}',
                  '--bias', f'{self.bias}',
                  '--extra_mem_dist', f'{self.extra_mem_dist}',
                  '--scaling_factor_start', f'{self.scaling_factor_start}',
                  '--scaling_factor_end', f'{self.scaling_factor_end}',
                  '--scaling_factor_step', f'{self.scaling_factor_step}',
                  '--cpu', f'{self.cpu}',
                  '--batch_size', f'{self.batch_size}']
        with open('run.out', 'w') as f:
            self.process = subprocess.Popen(['python','-u', '-m', 'memxterminator.radonfit.bin.membrane_subtract-main'] + params, stdout=f, stderr=subprocess.STDOUT)
        print("Membrane Subtraction started with PID:", self.process.pid)
        with open(self.PID_FILE, 'w') as f:
            f.write(str(self.process.pid))
    def kill_process(self):
        if self.process:
            self.process.terminate()
            print(f"Process PID {self.process.pid} terminated")
            self.process = None
            if os.path.exists(self.PID_FILE):
                os.remove(self.PID_FILE)
            self.timer.stop()
    def update_log(self):
        # 读取日志文件内容
        try:
            with open('run.out', 'r') as f:
                f.seek(self.last_read_position)  # 跳转到上次读取的位置
                new_content = f.read()  # 读取新内容
                self.last_read_position = f.tell()  # 更新读取的位置
            if new_content:
                self.textBrowser_log.append(new_content)
        except FileNotFoundError:
            self.textBrowser_log.append("Error: 'run.out' file not found.")
    def check_running_process(self):
        if os.path.exists(self.PID_FILE):
            with open(self.PID_FILE, 'r') as f:
                pid = int(f.read().strip())
            try:
                os.kill(pid, 0)  # Check if process is running
                self.process = pid
                print(f"Process with PID {pid} is still running!")
                # Here, you can also update the GUI to show that the process is running
            except OSError:
                print(f"Process with PID {pid} is not running.")
                os.remove(self.PID_FILE)


class MicrographMembraneSubtraction_Radon_App(QtWidgets.QDialog, Ui_MicrographMembraneSubtraction_Radonfit):
    PID_FILE = 'process.pid'
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.particles_selected_starfile_browse_pushButton.clicked.connect(self.particle_browse)
        self.process = None

        self.cpus = self.cpus_lineEdit.text()
        self.batch_size = self.batch_size_lineEdit.text()
        self.particle = None
        self.check_running_process()
        self.launch_pushButton.clicked.connect(self.start_process)
        self.kill_pushButton.clicked.connect(self.kill_process)
        
        with open('run.out', 'a') as f:
            f.write(f"Radonfit Micrograph Membrane Subtraction ready.\n")
        self.last_read_position = 0
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_log)
        self.timer.start(1000)
    
    def particle_browse(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open STAR File", "", "STAR Files (*.star)")
        if filepath:
            self.particles_selected_starfile_lineEdit.setText(filepath)
            self.particle = filepath
    
    def start_process(self):
        params = ['--particles_selected_filename', f'{self.particle}',
            '--cpu', f'{int(self.cpus)}',
            '--batch_size', f'{int(self.batch_size)}']
        with open('run.out', 'w') as f:
            self.process = subprocess.Popen(['python','-u', '-m', 'memxterminator.radonfit.bin.micrograph_mem_subtract'] + params, stdout=f, stderr=subprocess.STDOUT)
        print("Radon Micrograph Membrane Subtraction started with PID:", self.process.pid)
        with open(self.PID_FILE, 'w') as f:
            f.write(str(self.process.pid))
    def kill_process(self):
        if self.process:
            self.process.terminate()
            print(f"Process PID {self.process.pid} terminated")
            self.process = None
            if os.path.exists(self.PID_FILE):
                os.remove(self.PID_FILE)
            self.timer.stop()
    def update_log(self):
        # 读取日志文件内容
        try:
            with open('run.out', 'r') as f:
                f.seek(self.last_read_position)  # 跳转到上次读取的位置
                new_content = f.read()  # 读取新内容
                self.last_read_position = f.tell()  # 更新读取的位置
            if new_content:
                self.LOG_textBrowser.append(new_content)
        except FileNotFoundError:
            self.LOG_textBrowser.append("Error: 'run.out' file not found.")
    def check_running_process(self):
        if os.path.exists(self.PID_FILE):
            with open(self.PID_FILE, 'r') as f:
                pid = int(f.read().strip())
            try:
                os.kill(pid, 0)  # Check if process is running
                self.process = pid
                print(f"Process with PID {pid} is still running!")
                # Here, you can also update the GUI to show that the process is running
            except OSError:
                print(f"Process with PID {pid} is not running.")
                os.remove(self.PID_FILE)