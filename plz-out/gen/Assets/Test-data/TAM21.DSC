#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	NODATA
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	1024
XMIN	3488.700000
XWID	4.995117
*
*	Documentational Text:
*
TITL	'TAM 0.5mM 21%'
IRNAM	'Intensity'
XNAM	'Field'
IRUNI	''
XUNI	'G'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    irpe
DATE    01/24/13
TIME    14:33:19
CMNT    
SAMP    
SFOR    
STAG    L
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    NONE
AXS3    
A1CT    0.34912
A1SW    0.0005
MWFQ    9.812699e+09
MWPW    0.0009922
AVGS    1
RESO    TMS_0411
SPTP    0.01024
RCAG    10
RCHM    1
B0MA    5e-06
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    1024
RCTC    0.00512
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*
.DVC     acqStart, 1.0
.DVC     fieldCtrl, 1.0
AllegroMode        True
CenterField        3491.20 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
SetToSampleG       False
StaticFieldMon     3480.000 G
SweepDirection     Up
SweepWidth         5.0 G
.DVC     fieldSweep, 1.0
.DVC     freqCounter, 1.0
FrequencyMon       9.812699 GHz
QMonitBridge       On
.DVC     mwBridge, 1.0
AcqFineTuning      Each Slice Scan
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        100.0
Power              0.9922 mW
PowerAtten         23 dB
.DVC     recorder, 1.0
BaselineCorr       On
NbScansAcc         1
NbScansDone        1
NbScansToDo        1
ReplaceMode        Off
SmoothMode         Auto
SmoothPoints       5
.DVC     scanEnd, 1.0
.DVC     signalChannel, 1.0
AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvTime           10.24 ms
DModDetectSCT      First
DualDetect         OFF
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               10 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             0.050 G
ModFreq            100.00 kHz
ModPhase           0.0
Offset             0.0 %
QuadMode           False
Resolution         1024
Resonator          1
SctNorm            False
SctRevision        Allegro
SpuExtension       False
SweepTime          10.49 s
TimeConst          5.12 ms
TimeExp            False
TuneCaps           39
*
************************************************************
