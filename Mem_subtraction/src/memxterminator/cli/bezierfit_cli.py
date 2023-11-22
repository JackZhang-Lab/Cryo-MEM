from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QMainWindow, QFileDialog, QMessageBox, QApplication, QDialog
from memxterminator.GUI.MembraneAnalyzer_Bezierfit import Ui_MembraneAnalyzer_Bezierfit
from memxterminator.GUI.particle_membrane_subtraction_bezierfit import Ui_ParticleMembraneSubtraction_bezierfit
from memxterminator.GUI.MicrographMembraneSubtraction_Bezierfit import Ui_MicrographMembraneSubtraction_Bezierfit
import subprocess
import os
import json

class MembraneAnalyzer_Bezier_App(QtWidgets.QDialog, Ui_MembraneAnalyzer_Bezierfit):
    PID_FILE = 'process.pid'
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.particle_browse_pushButton.clicked.connect(self.particle_browse)
        self.template_browse_pushButton.clicked.connect(self.template_browse)
        self.output_pushButton.clicked.connect(self.output_browse)

        self.process = None
        self.template_filename = None
        self.particle_filename = None
        self.output_filename = None

        self.template_filename = self.template_lineEdit.text()
        self.particle_filename = self.particle_lineEdit.text()
        self.output_filename = self.output_lineEdit.text()

        self.degree = self.bezier_curve_degree_lineEdit.text()
        self.physical_membrane_dist = self.physical_membrane_distance_lineEdit.text()
        self.coarsefit_points_num = self.coarsefit_points_number_lineEdit.text()
        self.coarsefit_iter = self.coarsefit_iter_lineEdit.text()
        self.coarsefit_cpus = self.coarsefit_cpus_lineEdit.text()
        self.curve_penalty_thr = self.curve_penalty_thr_lineEdit.text()
        self.dithering_range = self.dithering_range_lineEdit.text()
        self.refine_iter = self.refine_iter_lineEdit.text()
        self.refine_cpus = self.refine_cpus_lineEdit.text()

        self.check_running_process()
        self.launch_pushButton.clicked.connect(self.start_process)
        self.kill_pushButton.clicked.connect(self.kill_process)
        
        with open('run.out', 'a') as f:
            f.write(f"Bezier Membrane Analysis ready.\n")
        self.last_read_position = 0
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_log)
        self.timer.start(1000)
    
    def template_browse(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open CS File", "", "CS Files (*.cs)")
        if filepath:
            self.template_lineEdit.setText(filepath)
            self.template_filename = filepath

    def particle_browse(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open CS File", "", "CS Files (*.cs)")
        if filepath:
            self.particle_lineEdit.setText(filepath)
            self.particle_filename = filepath

    def output_browse(self):
        filepath, _ = QFileDialog.getSaveFileName(self, "Save JSON File", "", "JSON Files (*.json)")
        if filepath:
            self.output_lineEdit.setText(filepath)
            self.output_filename = filepath
    
    def start_process(self):
        self.template_filename = self.template_lineEdit.text()
        self.particle_filename = self.particle_lineEdit.text()
        self.output_filename = self.output_lineEdit.text()

        self.degree = self.bezier_curve_degree_lineEdit.text()
        self.physical_membrane_dist = self.physical_membrane_distance_lineEdit.text()
        self.coarsefit_points_num = self.coarsefit_points_number_lineEdit.text()
        self.coarsefit_iter = self.coarsefit_iter_lineEdit.text()
        self.coarsefit_cpus = self.coarsefit_cpus_lineEdit.text()
        self.curve_penalty_thr = self.curve_penalty_thr_lineEdit.text()
        self.dithering_range = self.dithering_range_lineEdit.text()
        self.refine_iter = self.refine_iter_lineEdit.text()
        self.refine_cpus = self.refine_cpus_lineEdit.text()
        params = ['--template', f'{self.template_filename}',
            '--particle', f'{self.particle_filename}',
            '--output', f'{self.output_filename}',
            '--degree', f'{int(self.degree)}',
            '--physical_membrane_dist', f'{int(self.physical_membrane_dist)}',
            '--num_points', f'{int(self.coarsefit_points_num)}',
            '--coarsefit_iter', f'{int(self.coarsefit_iter)}',
            '--coarsefit_cpus', f'{int(self.coarsefit_cpus)}',
            '--cur_penalty_thr', f'{float(self.curve_penalty_thr)}',
            '--dithering_range', f'{int(self.dithering_range)}',
            '--refine_iter', f'{int(self.refine_iter)}',
            '--refine_cpus', f'{int(self.refine_cpus)}']
        
        with open('run.out', 'w') as f:
            self.process = subprocess.Popen(['python','-u', '-m', 'memxterminator.bezierfit.bin.mem_analyze_main'] + params, stdout=f, stderr=subprocess.STDOUT)
        print("Bezier Membrane Analysis started with PID:", self.process.pid)
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

class ParticleMembraneSubtraction_Bezier_App(QtWidgets.QDialog, Ui_ParticleMembraneSubtraction_bezierfit):
    PID_FILE = 'process.pid'
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.template_browse_pushButton.clicked.connect(self.template_browse)
        self.particle_browse_pushButton.clicked.connect(self.particle_browse)
        self.control_points_pushButton.clicked.connect(self.control_points_browse)
        self.process = None
        self.template = None
        self.particle = None
        self.control_points_filename = None

        self.physical_membrane_dist = self.physical_membrane_dist_lineEdit.text()
        self.points_step = self.points_step_lineEdit.text()

        self.check_running_process()
        self.launch_pushButton.clicked.connect(self.start_process)
        self.kill_pushButton.clicked.connect(self.kill_process)
        
        with open('run.out', 'a') as f:
            f.write(f"Bezier Particle Membrane Subtraction ready.\n")
        self.last_read_position = 0
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_log)
        self.timer.start(1000)
    
    def template_browse(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open CS File", "", "CS Files (*.cs)")
        if filepath:
            self.template_lineEdit.setText(filepath)
            self.template = filepath

    def particle_browse(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open CS File", "", "CS Files (*.cs)")
        if filepath:
            self.particle_lineEdit.setText(filepath)
            self.particle = filepath

    def control_points_browse(self):
        filepath, _ = QFileDialog.getOpenFileName(self, "Open JSON File", "", "JSON Files (*.json)")
        if filepath:
            self.controlpoints_lineEdit.setText(filepath)
            self.control_points_filename = filepath
    
    def start_process(self):
        self.template = self.template_lineEdit.text()
        self.particle = self.particle_lineEdit.text()
        self.control_points_filename = self.controlpoints_lineEdit.text()
        params = ['--template', f'{self.template}',
            '--particle', f'{self.particle}',
            '--control_points', f'{self.control_points_filename}',
            '--physical_membrane_dist', f'{int(self.physical_membrane_dist)}',
            '--points_step', f'{float(self.points_step)}']
        with open('run.out', 'w') as f:
            self.process = subprocess.Popen(['python','-u', '-m', 'memxterminator.bezierfit.bin.mem_subtract_main'] + params, stdout=f, stderr=subprocess.STDOUT)
        print("Bezier Particle Membrane Subtraction started with PID:", self.process.pid)
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


class MicrographMembraneSubtraction_Bezier_App(QtWidgets.QDialog, Ui_MicrographMembraneSubtraction_Bezierfit):
    PID_FILE = 'process.pid'
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)
        self.particle_browse_pushButton_3.clicked.connect(self.particle_browse)
        self.process = None

        self.particle = None
        self.cpus = self.cpus_lineEdit_3.text()
        self.batch_size = self.batch_size_lineEdit_3.text()

        self.check_running_process()
        self.launch_pushButton_3.clicked.connect(self.start_process)
        self.kill_pushButton_3.clicked.connect(self.kill_process)
        
        with open('run.out', 'a') as f:
            f.write(f"Bezier Micrograph Membrane Subtraction ready.\n")
        self.last_read_position = 0
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.update_log)
        self.timer.start(1000)
    
    def particle_browse(self):
        # filepath, _ = QFileDialog.getOpenFileName(self, "Open CS File", "", "CS Files (*.cs)")
        filepath, _ = QFileDialog.getOpenFileName(self, "Open STAR File", "", "STAR Files (*.star)")
        if filepath:
            self.particle_lineEdit_3.setText(filepath)
            self.particle = filepath
    
    def start_process(self):
        self.particle = self.particle_lineEdit_3.text()
        params = ['--particle', f'{self.particle}',
            '--cpus', f'{int(self.cpus)}',
            '--batch_size', f'{int(self.batch_size)}']
        with open('run.out', 'w') as f:
            self.process = subprocess.Popen(['python','-u', '-m', 'memxterminator.bezierfit.bin.micrograph_mem_subtract_main'] + params, stdout=f, stderr=subprocess.STDOUT)
        print("Bezier Micrograph Membrane Subtraction started with PID:", self.process.pid)
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