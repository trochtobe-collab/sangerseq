#!/usr/bin/env python3
"""
Sanger Sequencing Pipeline — GUI Launcher  (PyQt6 edition)
===========================================================
INSTALL (once):
    pip3 install PyQt6

USAGE:
    python3 sangerseq_gui.py
    (keep this file in the same folder as sangerseq-full2.py)
"""

import os, sys, json, subprocess, csv as _csv
from pathlib import Path

try:
    from PyQt6.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QTabWidget, QLabel, QPushButton, QLineEdit, QFileDialog,
        QTextEdit, QSlider, QScrollArea, QFrame, QSplitter,
        QListWidget, QListWidgetItem, QProgressBar, QSpinBox,
        QCheckBox, QComboBox, QTableWidget, QTableWidgetItem,
        QHeaderView, QAbstractItemView,
    )
    from PyQt6.QtCore import Qt, QThread, pyqtSignal
    from PyQt6.QtGui  import QFont, QColor, QTextCursor, QTextCharFormat
except ImportError:
    print("PyQt6 not found.\nInstall it with:  pip3 install PyQt6")
    sys.exit(1)

# ── Config (persists paths & params between sessions) ────────────────────────
CONFIG_PATH = Path.home() / ".sangerseq_gui.json"

def _load_config():
    try: return json.loads(CONFIG_PATH.read_text())
    except: return {}

def _save_config(data):
    try: CONFIG_PATH.write_text(json.dumps(data, indent=2))
    except: pass

# ── Parameter definitions ─────────────────────────────────────────────────────
PARAMS = [
    ("_SECTION_","TRIMMING","Section 1 — Read Quality Trimming",
     "Controls how low-quality bases are removed from the ENDS of each read.\nInterior bases are never trimmed."),
    ("Trim Window (bp)","TRIM_WINDOW",10,3,20,1,
     "How many bases the sliding window checks at once.\n\n↑ Higher → smoother, less sensitive to a single bad end base\n↓ Lower  → trims short noisy patches more aggressively\n\nRecommended: 10"),
    ("Trim Min Quality (Phred)","TRIM_MIN_QUAL",13,5,30,1,
     "The window average Phred must reach this before trimming stops.\n\nQ10=90%  Q13=95% ← rec  Q15=97%  Q20=99%\n\n↑ Higher → shorter reads, more inter-amplicon gaps\n↓ Lower  → longer reads, better amplicon overlap\n\nRecommended: 13"),
    ("Min Keep (bp)","MIN_KEEP_BP",500,100,1000,50,
     "Hard floor — trimming will NEVER reduce a read below this length.\nA WARNING is printed if this floor is enforced.\n\nRecommended: 500"),
    ("Min Usable (bp)","MIN_USABLE_BP",100,50,500,10,
     "Reads shorter than this after trimming are discarded entirely.\nOnly affects severely degraded samples.\n\nRecommended: 100"),
    ("_SECTION_","PILEUP","Section 2 — Pileup / Assembly Quality",
     "Controls how low-quality bases are treated when reads are stacked to build a consensus."),
    ("Min Deposit Quality (Phred)","MIN_DEPOSIT_QUAL",5,1,20,1,
     "Bases below this Phred are deposited as lowercase 'n' (nearly unreadable).\nThey still vote but LOSE to any confident uppercase call.\n\n↑ Higher → more bases demoted\n↓ Lower  → noisy bases keep their identity\n\nRecommended: 5"),
    ("Mask Quality (Phred)","MASK_QUAL",15,5,30,1,
     "Bases below this are hidden from the aligner during read anchoring.\nDoes NOT change the deposited base — only helps find the correct ref position.\n\n↑ Higher → cleaner anchor\n↓ Lower  → aligner sees noisier bases\n\nRecommended: 15"),
    ("_SECTION_","ALIGNMENT","Section 3 — Alignment Gap Penalties",
     "Controls how strongly the aligner resists opening gaps.\nHigh penalties prevent false frameshifts from sequencing artefacts."),
    ("Gap Open Penalty","ALIGN_GAP_OPEN",-20,-50,-5,1,
     "Penalty for opening a gap in the final contig-vs-reference alignment.\nMore negative = harder to open = fewer frameshifts.\n\n↑ Less negative → more gaps (phantom frameshifts risk)\n↓ More negative → almost no gaps\n\nRecommended: −20"),
    ("Gap Extend Penalty","ALIGN_GAP_EXTEND",-1,-10,-1,1,
     "Extra cost for each additional base in an open gap.\n\nRecommended: −1"),
    ("Anchor Ref Gap Open","ANCHOR_REF_GAP_OPEN",-5,-20,-1,1,
     "Penalty for a deletion in a read vs reference during Step 2.\nReal deletions are possible so only moderately penalised.\n\nRecommended: −5"),
    ("Anchor Query Gap Open","ANCHOR_QUERY_GAP_OPEN",-1000,-2000,-100,100,
     "Penalty for an insertion in a read vs reference during Step 2.\nSet very high to effectively forbid insertions (almost always artefacts).\nDo NOT change unless studying real insertions.\n\nRecommended: −1000"),
]
RECOMMENDED = {p[1]: p[2] for p in PARAMS if p[0] != "_SECTION_"}

PRESETS = {
    "Standard  (recommended)": dict(RECOMMENDED),
    "Strict  (high-quality data)": {
        "TRIM_WINDOW":10,"TRIM_MIN_QUAL":20,"MIN_KEEP_BP":500,"MIN_USABLE_BP":200,
        "MIN_DEPOSIT_QUAL":10,"MASK_QUAL":20,"ALIGN_GAP_OPEN":-25,"ALIGN_GAP_EXTEND":-1,
        "ANCHOR_REF_GAP_OPEN":-5,"ANCHOR_QUERY_GAP_OPEN":-1000,
    },
    "Lenient  (degraded samples)": {
        "TRIM_WINDOW":10,"TRIM_MIN_QUAL":10,"MIN_KEEP_BP":300,"MIN_USABLE_BP":100,
        "MIN_DEPOSIT_QUAL":3,"MASK_QUAL":10,"ALIGN_GAP_OPEN":-15,"ALIGN_GAP_EXTEND":-1,
        "ANCHOR_REF_GAP_OPEN":-5,"ANCHOR_QUERY_GAP_OPEN":-1000,
    },
}

C = {
    "bg":"#1e1e2e","panel":"#2a2a3e","entry":"#313244","border":"#45475a",
    "accent":"#7c6af7","teal":"#56cfb2","green":"#3ddc84",
    "orange":"#f5a623","red":"#ff5f5f","fg":"#cdd6f4","fg_dim":"#7f849c",
    "log_bg":"#11111b",
}

STYLE = f"""
QMainWindow,QWidget{{background:{C['bg']};color:{C['fg']};}}
QTabWidget::pane{{border:1px solid {C['border']};background:{C['bg']};}}
QTabBar::tab{{background:{C['panel']};color:{C['fg_dim']};padding:9px 22px;border:none;border-top-left-radius:6px;border-top-right-radius:6px;font-size:12px;}}
QTabBar::tab:selected{{background:{C['accent']};color:white;font-weight:bold;}}
QTabBar::tab:hover{{background:{C['entry']};color:{C['fg']};}}
QLineEdit{{background:{C['entry']};color:{C['fg']};border:1px solid {C['border']};border-radius:5px;padding:5px 8px;font-size:12px;}}
QLineEdit:focus{{border:1px solid {C['accent']};}}
QPushButton{{background:{C['entry']};color:{C['fg']};border:1px solid {C['border']};border-radius:6px;padding:6px 14px;font-size:12px;}}
QPushButton:hover{{background:{C['border']};}}
QPushButton:disabled{{background:{C['border']};color:{C['fg_dim']};}}
QScrollBar:vertical{{background:{C['panel']};width:8px;border-radius:4px;}}
QScrollBar::handle:vertical{{background:{C['border']};border-radius:4px;min-height:20px;}}
QScrollBar::add-line:vertical,QScrollBar::sub-line:vertical{{height:0px;}}
QScrollBar:horizontal{{background:{C['panel']};height:8px;border-radius:4px;}}
QScrollBar::handle:horizontal{{background:{C['border']};border-radius:4px;min-width:20px;}}
QScrollBar::add-line:horizontal,QScrollBar::sub-line:horizontal{{width:0px;}}
QListWidget{{background:{C['entry']};color:{C['fg']};border:1px solid {C['border']};border-radius:6px;font-family:Menlo,Consolas,monospace;font-size:11px;}}
QListWidget::item:selected{{background:{C['accent']};color:white;}}
QListWidget::item:hover{{background:{C['panel']};}}
QLabel{{color:{C['fg']};background:transparent;}}
QToolTip{{background:{C['panel']};color:{C['fg']};border:1px solid {C['border']};padding:8px;font-size:12px;}}
QProgressBar{{background:{C['entry']};border:none;border-radius:4px;}}
QProgressBar::chunk{{background:{C['accent']};border-radius:4px;}}
QSplitter::handle{{background:{C['border']};width:2px;height:2px;}}
QSpinBox{{background:{C['panel']};color:{C['fg']};border:1px solid {C['border']};border-radius:5px;padding:4px;font-size:12px;}}
QCheckBox{{color:{C['fg_dim']};font-size:11px;spacing:6px;}}
QCheckBox::indicator{{width:14px;height:14px;border-radius:3px;border:1px solid {C['border']};background:{C['entry']};}}
QCheckBox::indicator:checked{{background:{C['accent']};border-color:{C['accent']};}}
QComboBox{{background:{C['entry']};color:{C['fg']};border:1px solid {C['border']};border-radius:5px;padding:5px 8px;font-size:12px;}}
QComboBox::drop-down{{border:none;width:20px;}}
QComboBox QAbstractItemView{{background:{C['panel']};color:{C['fg']};border:1px solid {C['border']};selection-background-color:{C['accent']};}}
QTableWidget{{background:{C['entry']};color:{C['fg']};border:1px solid {C['border']};border-radius:6px;gridline-color:{C['border']};font-size:11px;}}
QTableWidget::item{{padding:4px 8px;}}
QTableWidget::item:selected{{background:{C['accent']};color:white;}}
QHeaderView::section{{background:{C['panel']};color:{C['accent']};border:none;border-bottom:1px solid {C['border']};padding:6px 8px;font-size:11px;font-weight:bold;}}
"""

# ── Worker thread ─────────────────────────────────────────────────────────────
class PipelineWorker(QThread):
    log_line = pyqtSignal(str)
    finished = pyqtSignal(bool)

    def __init__(self, raw_ab1, ref_a, ref_b, output_dir, params):
        super().__init__()
        self.raw_ab1,self.ref_a,self.ref_b = raw_ab1,ref_a,ref_b
        self.output_dir,self.params = output_dir,params
        self._proc,self._stop = None,False

    def run(self):
        gui_dir = Path(__file__).parent
        pipeline_py = gui_dir/"sangerseq-full2.py"
        if not pipeline_py.exists():
            self.log_line.emit(f"ERROR: Cannot find sangerseq-full2.py in {gui_dir}")
            self.finished.emit(False); return
        runner = gui_dir/"_gui_runner_tmp.py"
        p,od = self.params,self.output_dir
        neg = "{'ALIGN_GAP_OPEN','ALIGN_GAP_EXTEND','ANCHOR_REF_GAP_OPEN','ANCHOR_QUERY_GAP_OPEN'}"
        code = f"""
import os,sys,shutil; sys.path.insert(0,r'{gui_dir}')
import importlib.util
spec=importlib.util.spec_from_file_location("sangerseq",r'{pipeline_py}')
mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
for name,val in {p!r}.items():
    neg={neg}
    final=-abs(int(val)) if name in neg else int(val)
    setattr(mod,name,final); print(f'[GUI] {{name}} = {{final}}')
from Bio import SeqIO
RAW=r'{self.raw_ab1}'; RFA=r'{self.ref_a}'; RFB=r'{self.ref_b}'; OUT=r'{od}'
F1=os.path.join(OUT,'01_amplicon_reads'); F2=os.path.join(OUT,'02_assembled_contigs')
F3=os.path.join(OUT,'03_analysis_reports'); F4=os.path.join(OUT,'04_extracted_F_protein')
CSV=os.path.join(F3,'master_mutation_summary.csv')
print('='*60+'\\n  Sanger Sequencing Pipeline\\n'+'='*60)
try:
    rA=SeqIO.read(RFA,'fasta'); rB=SeqIO.read(RFB,'fasta')
    print(f'INFO: Ref A: {{rA.id}} ({{len(rA.seq)}} bp)'); print(f'INFO: Ref B: {{rB.id}} ({{len(rB.seq)}} bp)')
except Exception as e: print(f'ERROR: {{e}}'); sys.exit(1)
refs={{'RSV_A':rA,'RSV_B':rB}}
os.makedirs(OUT,exist_ok=True)
for f in [F3,F4,F1,F2]:
    if os.path.exists(f): shutil.rmtree(f)
    os.makedirs(f,exist_ok=True)
ms=[]
print('\\n[STEP 1] Reading chromatograms...')
amps=mod.step1_process_ab1_folder(RAW,F1)
for strain,samps in amps.items():
    if not samps: print(f'INFO: No {{strain}} — skipping.'); continue
    ref=refs[strain]
    print(f'\\n[STEP 2] Assembling {{strain}} ({{len(samps)}} samples)...')
    cp=mod.step2_assemble_contigs(samps,os.path.join(F2,strain),ref)
    print(f'\\n[STEPS 3-5] Calling mutations — {{strain}}...')
    mod.process_batch(ref,cp,os.path.join(F3,strain),os.path.join(F4,strain),strain,ms)
if ms:
    print('\\n[OUTPUT] Writing master CSV...'); mod.write_master_csv(ms,CSV)
print('\\n[OUTPUT] Writing summary...')
mod.write_pipeline_summary(summary_path=os.path.join(OUT,'pipeline_summary.txt'),master_summary=ms,ref_a=rA,ref_b=rB)
print('\\n[DONE] Pipeline complete.')
"""
        runner.write_text(code)
        try:
            self._proc=subprocess.Popen([sys.executable,str(runner)],
                stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
            for line in self._proc.stdout:
                if self._stop: break
                self.log_line.emit(line.rstrip())
            self._proc.wait()
            self.finished.emit((self._proc.returncode==0) and not self._stop)
        except Exception as ex:
            self.log_line.emit(f"ERROR: {ex}"); self.finished.emit(False)
        finally:
            try: runner.unlink()
            except: pass

    def cancel(self):
        self._stop=True
        if self._proc: self._proc.terminate()

# ── Widget helpers ────────────────────────────────────────────────────────────
def _card():
    f=QFrame()
    f.setStyleSheet(f"QFrame{{background:{C['panel']};border-radius:10px;border:1px solid {C['border']};}}")
    return f

def _accent_btn(text,color,text_color="white"):
    b=QPushButton(text)
    b.setStyleSheet(f"QPushButton{{background:{color};color:{text_color};border:none;border-radius:6px;padding:8px 18px;font-size:13px;font-weight:bold;}}QPushButton:hover{{background:{color}cc;}}QPushButton:disabled{{background:{C['border']};color:{C['fg_dim']};}}")
    return b

def _dim_btn(text):
    b=QPushButton(text)
    b.setStyleSheet(f"QPushButton{{background:{C['entry']};color:{C['fg_dim']};border:1px solid {C['border']};border-radius:5px;padding:5px 12px;font-size:11px;}}QPushButton:hover{{color:{C['fg']};border-color:{C['accent']};}}")
    return b

def _sec_lbl(text):
    l=QLabel(text)
    l.setStyleSheet(f"color:{C['accent']};font-size:11px;font-weight:bold;letter-spacing:1px;padding:4px 0;")
    return l

def _mono():
    t=QTextEdit(); t.setReadOnly(True)
    t.setFont(QFont("Menlo",11))
    t.setStyleSheet(f"QTextEdit{{background:{C['log_bg']};color:{C['fg']};border:none;border-radius:6px;padding:8px;}}")
    return t

def _err_lbl():
    l=QLabel(""); l.setStyleSheet(f"color:{C['red']};font-size:11px;padding-left:152px;"); l.setVisible(False); return l

# ── Main window ───────────────────────────────────────────────────────────────
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("🧬  Sanger Sequencing Pipeline")
        self.resize(1200,860); self.setMinimumSize(960,680)
        self.setStyleSheet(STYLE)
        self._cfg=_load_config(); self._worker=None; self._running=False
        self._sliders={}; self._spins={}; self._report_paths=[]; self._current_report=None
        self._all_report_items=[]; self._step=0

        central=QWidget(); self.setCentralWidget(central)
        vl=QVBoxLayout(central); vl.setContentsMargins(0,0,0,0); vl.setSpacing(0)
        vl.addWidget(self._header())
        self.tabs=QTabWidget(); self.tabs.setDocumentMode(True)
        self.tabs.addTab(self._run_tab(),    "  ▶  Run  ")
        self.tabs.addTab(self._params_tab(), "  ⚙  Parameters  ")
        self.tabs.addTab(self._results_tab(),"  📋  Results  ")
        vl.addWidget(self.tabs)

    def _header(self):
        h=QWidget(); h.setFixedHeight(56)
        h.setStyleSheet(f"background:{C['panel']};border-bottom:1px solid {C['border']};")
        hl=QHBoxLayout(h); hl.setContentsMargins(20,0,20,0)
        t=QLabel("🧬  Sanger Sequencing Pipeline"); t.setStyleSheet(f"color:{C['fg']};font-size:18px;font-weight:bold;"); hl.addWidget(t)
        s=QLabel("RSV-A / RSV-B  ·  Full Automation  ·  Steps 1–5"); s.setStyleSheet(f"color:{C['fg_dim']};font-size:12px;"); hl.addWidget(s); hl.addStretch(); return h

    # ── Run tab ───────────────────────────────────────────────────────────────
    def _run_tab(self):
        w=QWidget(); lay=QVBoxLayout(w); lay.setContentsMargins(14,10,14,10); lay.setSpacing(8)
        # Paths card
        pc=_card(); pl=QVBoxLayout(pc); pl.setContentsMargins(14,10,14,12); pl.setSpacing(2)
        pl.addWidget(_sec_lbl("INPUT / OUTPUT PATHS"))
        self._paths={}; self._path_errs={}
        for label,key,default,is_dir in [
            ("Raw .ab1 folder","raw_ab1","input/Raw_Data",True),
            ("Reference RSV-A","ref_a","input/references/ref_RSV_A.fasta",False),
            ("Reference RSV-B","ref_b","input/references/ref_RSV_B.fasta",False),
            ("Output folder","output_dir","output",True),
        ]:
            saved=self._cfg.get(f"path_{key}",default)
            row=QHBoxLayout(); lbl=QLabel(label); lbl.setFixedWidth(145)
            lbl.setStyleSheet(f"color:{C['fg_dim']};font-size:12px;"); row.addWidget(lbl)
            edit=QLineEdit(saved); self._paths[key]=edit
            edit.textChanged.connect(lambda _,k=key: self._validate_path(k)); row.addWidget(edit)
            btn=_dim_btn("Browse"); btn.setFixedWidth(72)
            btn.clicked.connect(lambda _,k=key,d=is_dir: self._browse(k,d)); row.addWidget(btn)
            pl.addLayout(row)
            err=_err_lbl(); self._path_errs[key]=err; pl.addWidget(err)
        open_out=_dim_btn("📂  Open output folder in Finder"); open_out.clicked.connect(self._open_output_folder); pl.addWidget(open_out)
        lay.addWidget(pc)
        # Controls
        ctrl=QHBoxLayout()
        self.run_btn=_accent_btn("▶   Run Pipeline",C["green"],"#1e1e2e"); self.run_btn.setFixedSize(190,42); self.run_btn.clicked.connect(self._toggle_run); ctrl.addWidget(self.run_btn)
        self.status_lbl=QLabel("Ready"); self.status_lbl.setStyleSheet(f"color:{C['fg_dim']};font-size:12px;padding-left:10px;"); ctrl.addWidget(self.status_lbl); ctrl.addStretch()
        pcol=QVBoxLayout(); pcol.setSpacing(2)
        self.progress=QProgressBar(); self.progress.setFixedSize(300,8); self.progress.setTextVisible(False); self.progress.setValue(0); self.progress.setRange(0,4); pcol.addWidget(self.progress)
        self.step_lbl=QLabel(""); self.step_lbl.setStyleSheet(f"color:{C['fg_dim']};font-size:10px;"); pcol.addWidget(self.step_lbl); ctrl.addLayout(pcol)
        lay.addLayout(ctrl)
        # Log card
        lc=_card(); ll=QVBoxLayout(lc); ll.setContentsMargins(12,8,12,12); ll.setSpacing(4)
        lt=QHBoxLayout(); lt.addWidget(_sec_lbl("LIVE LOG")); lt.addStretch()
        self.autoscroll_chk=QCheckBox("Auto-scroll"); self.autoscroll_chk.setChecked(True); lt.addWidget(self.autoscroll_chk)
        clr=_dim_btn("Clear"); clr.clicked.connect(self._clear_log); lt.addWidget(clr); ll.addLayout(lt)
        self.log_edit=_mono(); ll.addWidget(self.log_edit); lay.addWidget(lc,stretch=1)
        # Summary card
        self.summary_card=_card(); self.summary_card.setVisible(False)
        sl=QVBoxLayout(self.summary_card); sl.setContentsMargins(14,10,14,12); sl.setSpacing(4)
        sl.addWidget(_sec_lbl("LAST RUN SUMMARY"))
        self.summary_text=QLabel(""); self.summary_text.setStyleSheet(f"color:{C['fg']};font-size:12px;"); self.summary_text.setWordWrap(True); sl.addWidget(self.summary_text)
        lay.addWidget(self.summary_card)
        return w

    # ── Params tab ────────────────────────────────────────────────────────────
    def _params_tab(self):
        outer=QWidget(); ol=QVBoxLayout(outer); ol.setContentsMargins(14,10,14,10); ol.setSpacing(6)
        top=QHBoxLayout()
        top.addWidget(QLabel("Preset:"))
        self.preset_combo=QComboBox(); self.preset_combo.addItem("— select a preset —")
        for name in PRESETS: self.preset_combo.addItem(name)
        self.preset_combo.setFixedWidth(260); self.preset_combo.currentTextChanged.connect(self._apply_preset); top.addWidget(self.preset_combo); top.addStretch()
        info=QLabel("Hover over a parameter name for details."); info.setStyleSheet(f"color:{C['fg_dim']};font-size:12px;"); top.addWidget(info)
        rst=_dim_btn("Reset all to recommended"); rst.clicked.connect(self._reset_params); top.addWidget(rst); ol.addLayout(top)
        sa=QScrollArea(); sa.setWidgetResizable(True); sa.setStyleSheet(f"QScrollArea{{background:{C['bg']};border:none;}}")
        content=QWidget(); content.setStyleSheet(f"background:{C['bg']};"); cl=QVBoxLayout(content); cl.setContentsMargins(0,0,8,0); cl.setSpacing(4)
        for entry in PARAMS:
            if entry[0]=="_SECTION_":
                sec=_card(); sec.setStyleSheet(f"QFrame{{background:{C['panel']};border-radius:8px;border:1px solid {C['border']};}}"); sl=QVBoxLayout(sec); sl.setContentsMargins(14,8,14,8); sl.setSpacing(1)
                tl=QLabel(entry[2].upper()); tl.setStyleSheet(f"color:{C['accent']};font-size:12px;font-weight:bold;"); sl.addWidget(tl)
                dl=QLabel(entry[3]); dl.setStyleSheet(f"color:{C['fg_dim']};font-size:11px;"); dl.setWordWrap(True); sl.addWidget(dl)
                cl.addWidget(sec); cl.addSpacing(2); continue
            label,name,default,vmin,vmax,step,tooltip=entry
            rw=QWidget(); rw.setStyleSheet(f"QWidget{{background:{C['entry']};border-radius:6px;}}"); rl=QHBoxLayout(rw); rl.setContentsMargins(12,6,12,6); rl.setSpacing(10)
            lbl=QLabel(label); lbl.setFixedWidth(210); lbl.setStyleSheet(f"color:{C['fg']};font-size:12px;background:transparent;"); lbl.setToolTip(tooltip); lbl.setCursor(Qt.CursorShape.WhatsThisCursor); rl.addWidget(lbl)
            slider=QSlider(Qt.Orientation.Horizontal); slider.setMinimum(vmin); slider.setMaximum(vmax); slider.setValue(self._cfg.get(f"param_{name}",default)); slider.setSingleStep(step)
            slider.setStyleSheet(f"QSlider::groove:horizontal{{background:{C['border']};height:4px;border-radius:2px;}}QSlider::handle:horizontal{{background:{C['accent']};width:14px;height:14px;margin:-5px 0;border-radius:7px;}}QSlider::sub-page:horizontal{{background:{C['accent']};border-radius:2px;}}")
            self._sliders[name]=slider; rl.addWidget(slider,stretch=1)
            spin=QSpinBox(); spin.setMinimum(vmin); spin.setMaximum(vmax); spin.setValue(self._cfg.get(f"param_{name}",default)); spin.setSingleStep(step); spin.setFixedWidth(82); self._spins[name]=spin; rl.addWidget(spin)
            slider.valueChanged.connect(lambda v,s=spin:   s.setValue(v))
            spin.valueChanged.connect(  lambda v,s=slider: s.setValue(v))
            spin.valueChanged.connect(  lambda v,n=name:   self._save_param(n,v))
            rec=QLabel(f"rec: {default}"); rec.setFixedWidth(76); rec.setStyleSheet(f"color:{C['fg_dim']};font-size:10px;background:transparent;"); rl.addWidget(rec); cl.addWidget(rw)
        cl.addStretch(); sa.setWidget(content); ol.addWidget(sa,stretch=1); return outer

    # ── Results tab ───────────────────────────────────────────────────────────
    def _results_tab(self):
        w=QWidget(); lay=QVBoxLayout(w); lay.setContentsMargins(14,10,14,10); lay.setSpacing(8)
        top=QHBoxLayout()
        info=QLabel("Click a report to view it.  CSV files are shown as a table."); info.setStyleSheet(f"color:{C['fg_dim']};font-size:12px;"); top.addWidget(info); top.addStretch()
        filter_lbl=QLabel("Filter:"); filter_lbl.setStyleSheet(f"color:{C['fg_dim']};font-size:12px;"); top.addWidget(filter_lbl)
        self.filter_edit=QLineEdit(); self.filter_edit.setFixedWidth(160); self.filter_edit.setPlaceholderText("sample name…"); self.filter_edit.textChanged.connect(self._filter_results); top.addWidget(self.filter_edit)
        rb=_dim_btn("↻  Refresh"); rb.clicked.connect(self._refresh_results); top.addWidget(rb); lay.addLayout(top)
        spl=QSplitter(Qt.Orientation.Horizontal)
        lc=_card(); ll2=QVBoxLayout(lc); ll2.setContentsMargins(10,10,10,10); ll2.setSpacing(4)
        ll2.addWidget(_sec_lbl("REPORT FILES"))
        self.file_list=QListWidget(); ll2.addWidget(self.file_list); self.file_list.currentRowChanged.connect(self._on_file_select); spl.addWidget(lc)
        vc=_card(); vl2=QVBoxLayout(vc); vl2.setContentsMargins(10,10,10,10); vl2.setSpacing(4)
        vt=QHBoxLayout(); self.viewer_title=QLabel("Select a report →"); self.viewer_title.setStyleSheet(f"color:{C['fg_dim']};font-size:12px;"); vt.addWidget(self.viewer_title); vt.addStretch()
        fb=_dim_btn("Open in Finder"); fb.clicked.connect(self._open_in_finder); vt.addWidget(fb)
        fob=_dim_btn("📂  Open folder"); fob.clicked.connect(self._open_output_folder); vt.addWidget(fob); vl2.addLayout(vt)
        self.viewer=_mono()
        self.table_view=QTableWidget()
        self.table_view.setFont(QFont("Menlo",11))
        self.table_view.setStyleSheet(f"QTableWidget{{background:{C['log_bg']};border-radius:6px;}}")
        self.table_view.setAlternatingRowColors(True)
        self.table_view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table_view.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table_view.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table_view.setVisible(False)
        vl2.addWidget(self.viewer,stretch=1); vl2.addWidget(self.table_view,stretch=1); spl.addWidget(vc); spl.setSizes([280,860])
        lay.addWidget(spl,stretch=1); return w

    # ── Logic ─────────────────────────────────────────────────────────────────
    def _browse(self,key,is_dir):
        start=self._paths[key].text() or str(Path.home())
        p=(QFileDialog.getExistingDirectory(self,"Select folder",start) if is_dir
           else QFileDialog.getOpenFileName(self,"Select file",start,"FASTA (*.fasta *.fa);;All (*)")[0])
        if p:
            self._paths[key].setText(p); self._cfg[f"path_{key}"]=p; _save_config(self._cfg)

    def _validate_path(self,key):
        text=self._paths[key].text().strip(); err=self._path_errs[key]
        if not text: err.setText("  ⚠  Path cannot be empty"); err.setVisible(True)
        elif not os.path.exists(text): err.setText(f"  ⚠  Not found: {text}"); err.setVisible(True)
        else: err.setVisible(False)

    def _all_paths_valid(self):
        return all(os.path.exists(self._paths[k].text().strip()) for k in self._paths)

    def _reset_params(self):
        for name,slider in self._sliders.items(): slider.setValue(int(RECOMMENDED[name]))

    def _apply_preset(self,text):
        if text not in PRESETS: return
        for name,val in PRESETS[text].items():
            if name in self._sliders: self._sliders[name].setValue(int(val))

    def _save_param(self,name,val):
        self._cfg[f"param_{name}"]=val; _save_config(self._cfg)

    def _toggle_run(self):
        if self._running:
            if self._worker: self._worker.cancel()
            self._done(cancelled=True)
        else: self._start()

    def _start(self):
        for key in self._paths: self._validate_path(key)
        if not self._all_paths_valid():
            self._set_status("Fix the path errors above before running.",C["red"]); return
        self._running=True; self._step=0
        self.run_btn.setText("⏹   Stop")
        self.run_btn.setStyleSheet(self.run_btn.styleSheet().replace(C["green"],C["red"]).replace("#1e1e2e","white"))
        self._set_status("Running…",C["orange"]); self.progress.setRange(0,4); self.progress.setValue(0)
        self.step_lbl.setText("Initialising…"); self.summary_card.setVisible(False)
        params={n:s.value() for n,s in self._spins.items()}
        self._worker=PipelineWorker(self._paths["raw_ab1"].text(),self._paths["ref_a"].text(),self._paths["ref_b"].text(),self._paths["output_dir"].text(),params)
        self._worker.log_line.connect(self._append_log); self._worker.finished.connect(self._on_finished); self._worker.start()

    def _on_finished(self,ok):
        if ok:
            self._append_log("\n✅  Pipeline finished successfully!")
            self._refresh_results(); self._load_run_summary(); self.tabs.setCurrentIndex(2)
        else:
            self._append_log("\n❌  Pipeline finished with errors — check log above.")
        self._done(error=not ok)

    def _done(self,cancelled=False,error=False):
        self._running=False
        self.run_btn.setText("▶   Run Pipeline")
        self.run_btn.setStyleSheet(self.run_btn.styleSheet().replace(C["red"],C["green"]).replace("white","#1e1e2e"))
        self.progress.setRange(0,4); self.progress.setValue(0 if (error or cancelled) else 4)
        self.step_lbl.setText("Stopped" if cancelled else ("Error" if error else "Complete ✓"))
        self._set_status("Stopped." if cancelled else ("Finished with errors." if error else "Done ✓"),
                         C["fg_dim"] if cancelled else (C["red"] if error else C["green"]))

    def _set_status(self,text,color=None):
        self.status_lbl.setText(text)
        self.status_lbl.setStyleSheet(f"color:{color or C['fg_dim']};font-size:12px;padding-left:10px;")

    def _append_log(self,line):
        lu=line.upper()
        if "[STEP 1]" in lu:   self.progress.setValue(1); self.step_lbl.setText("Step 1 / 4  —  Reading chromatograms")
        elif "[STEP 2]" in lu: self.progress.setValue(2); self.step_lbl.setText("Step 2 / 4  —  Assembling contigs")
        elif "[STEPS 3" in lu: self.progress.setValue(3); self.step_lbl.setText("Step 3–5 / 4  —  Calling mutations")
        elif "[OUTPUT]" in lu: self.progress.setValue(4); self.step_lbl.setText("Writing output files…")
        cur=self.log_edit.textCursor(); cur.movePosition(QTextCursor.MoveOperation.End)
        fmt=QTextCharFormat()
        if "ERROR" in lu or "CRITICAL" in lu:     fmt.setForeground(QColor(C["red"]))
        elif "WARNING" in lu:                      fmt.setForeground(QColor(C["orange"]))
        elif lu.startswith("[STEP") or lu.startswith("[OUTPUT") or lu.startswith("[DONE"):
            fmt.setForeground(QColor(C["teal"])); fmt.setFontWeight(700)
        elif lu.startswith("INFO") or "✅" in line: fmt.setForeground(QColor(C["green"]))
        elif "[GUI]" in line:                      fmt.setForeground(QColor(C["fg_dim"]))
        else:                                      fmt.setForeground(QColor(C["fg"]))
        cur.insertText(line+"\n",fmt)
        if self.autoscroll_chk.isChecked():
            self.log_edit.setTextCursor(cur); self.log_edit.ensureCursorVisible()

    def _clear_log(self): self.log_edit.clear()

    def _open_output_folder(self):
        od=self._paths["output_dir"].text()
        if os.path.isdir(od): subprocess.run(["open",od])

    def _load_run_summary(self):
        od=self._paths["output_dir"].text()
        csv_path=os.path.join(od,"03_analysis_reports","master_mutation_summary.csv")
        if not os.path.exists(csv_path): return
        try:
            with open(csv_path,newline="",encoding="utf-8") as f: rows=list(_csv.DictReader(f))
        except: return
        if not rows: return
        n=len(rows); na=sum(1 for r in rows if r.get("Strain")=="RSV_A"); nb=n-na
        muts=sum(int(r.get("Total_Mutations",0)) for r in rows)
        warns=sum(1 for r in rows if int(r.get("Residual_Stops",0))>0)
        warn_str=(f"<span style='color:{C['red']}'><b>⛔  {warns} sample(s) need manual review</b></span>" if warns
                  else f"<span style='color:{C['green']}'>✓  No unresolved stop codons</span>")
        self.summary_text.setText(
            f"<b>Samples:</b> {n}  (RSV-A: {na}  |  RSV-B: {nb})  &nbsp;&nbsp;"
            f"<b>Total AA mutations:</b> {muts}  &nbsp;&nbsp;  {warn_str}")
        self.summary_card.setVisible(True)

    def _refresh_results(self):
        od=self._paths["output_dir"].text(); self._all_report_items=[]
        if not os.path.isdir(od): return
        for root,dirs,files in os.walk(od):
            dirs[:] = [d for d in sorted(dirs) if not d.startswith(".")]
            for f in sorted(files):
                if f.endswith((".txt",".csv")):
                    full=os.path.join(root,f); rel=os.path.relpath(full,od)
                    self._all_report_items.append((rel,full))
        self._populate_file_list(self._all_report_items)

    def _populate_file_list(self,items):
        self.file_list.clear(); self._report_paths=[]
        for rel,full in items:
            self._report_paths.append(full)
            item=QListWidgetItem(rel)
            if rel.endswith(".csv"): item.setForeground(QColor(C["teal"]))
            self.file_list.addItem(item)

    def _filter_results(self,text):
        text=text.lower()
        self._populate_file_list([(r,f) for r,f in self._all_report_items if text in r.lower()])

    def _on_file_select(self,row):
        if row<0 or row>=len(self._report_paths): return
        path=self._report_paths[row]; self._current_report=path
        self.viewer_title.setText(os.path.basename(path))
        self.viewer_title.setStyleSheet(f"color:{C['fg']};font-size:12px;")
        if path.endswith(".csv"):
            self._show_csv_table(path); return
        self.table_view.setVisible(False); self.viewer.setVisible(True)
        try: text=Path(path).read_text(encoding="utf-8",errors="replace")
        except Exception as e: text=f"Could not read: {e}"
        self.viewer.clear()
        cur=self.viewer.textCursor(); cur.movePosition(QTextCursor.MoveOperation.End)
        for line in text.splitlines(keepends=True):
            s=line.strip(); fmt=QTextCharFormat()
            if s.startswith("═") or s.startswith("─"):
                fmt.setForeground(QColor(C["teal"])); fmt.setFontWeight(700)
            elif s.startswith("SAMPLE") or s.startswith("STRAIN") or s.startswith("REFERENCE"):
                fmt.setForeground(QColor(C["accent"])); fmt.setFontWeight(700)
            elif "ERROR" in s.upper() or s.startswith("!") or "STOP CODON" in s.upper() or "⛔" in s:
                fmt.setForeground(QColor(C["red"]))
            elif s.startswith("?") or "UNCERTAIN" in s.upper() or "IUPAC" in s.upper() or "⚠" in s:
                fmt.setForeground(QColor(C["orange"]))
            elif s.startswith("-") and "AA Position" in s:
                fmt.setForeground(QColor(C["fg"])); fmt.setFontWeight(600)
            elif "✓" in s or "►" in s:
                fmt.setForeground(QColor(C["green"]))
            elif s.startswith("[") and "]" in s:
                fmt.setForeground(QColor(C["accent"])); fmt.setFontWeight(700)
            else: fmt.setForeground(QColor(C["fg"]))
            cur.insertText(line,fmt)
        self.viewer.setTextCursor(cur)
        self.viewer.moveCursor(QTextCursor.MoveOperation.Start)

    def _show_csv_table(self,path):
        self.viewer.setVisible(False); self.table_view.setVisible(True)
        try:
            with open(path,newline="",encoding="utf-8") as f:
                reader=_csv.DictReader(f); rows=list(reader); cols=reader.fieldnames or []
        except Exception as e:
            self.viewer.setVisible(True); self.table_view.setVisible(False)
            self.viewer.setPlainText(f"Could not parse CSV: {e}"); return
        self.table_view.setRowCount(len(rows)); self.table_view.setColumnCount(len(cols))
        self.table_view.setHorizontalHeaderLabels(cols)
        for ri,row in enumerate(rows):
            for ci,col in enumerate(cols):
                val=row.get(col,""); item=QTableWidgetItem(val)
                item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                cl=col.lower()
                if cl=="total_mutations" and val.isdigit():
                    item.setForeground(QColor(C["orange"] if int(val)>5 else C["fg"]))
                elif cl=="residual_stops" and val!="0":
                    item.setForeground(QColor(C["red"])); item.setFont(QFont("Menlo",11,QFont.Weight.Bold))
                elif "identity" in cl:
                    try:
                        pct=float(val)
                        item.setForeground(QColor(C["orange"] if pct<95 else (C["green"] if pct>=99 else C["fg"])))
                    except: pass
                self.table_view.setItem(ri,ci,item)
        self.table_view.resizeColumnsToContents()

    def _open_in_finder(self):
        if self._current_report and os.path.exists(self._current_report):
            subprocess.run(["open","-R",self._current_report])

if __name__=="__main__":
    app=QApplication(sys.argv); app.setStyle("Fusion")
    win=MainWindow(); win.show(); sys.exit(app.exec())