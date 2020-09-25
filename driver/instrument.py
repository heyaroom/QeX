from retrying import retry
from .util import name_to_alpha
from measurement_tool.units import MHz, GHz, ns, us, dB

def set_qubit(instrument, qubit_notes, qubit_information):
    """Setup of the instruments used in the experiment
    Args:
        instrument (TimeDomainInstrumentManager): Instruments used in the experiment
        qubit_notes (AttributionDict): Calibration notes for single-qubit control
        qubit_information (list): List of the qubit names
    """
    for qubit_name in qubit_information:
        note                                              = qubit_notes[qubit_name]
        instrument.pma.config[qubit_name+"_readout"].freq = note.cavity_readout_window_frequency
        instrument.pma.config[qubit_name+"_qubit"].freq   = note.qubit_dressed_frequency_jazz
        instrument.mwb.config[qubit_name].ratten          = note.cavity_readout_attenuation
        instrument.mwb.config[qubit_name].qatten          = note.pi_pulse_qubit_pump_attenuation
        instrument.mwb.config[qubit_name].rxatten         = 0*dB
        instrument.mwb.config[qubit_name].rgain           = note.cavity_readout_gain
        instrument.mwb.config[qubit_name].qgain           = note.pi_pulse_qubit_pump_gain
        instrument.mwb.config[qubit_name].c1gain          = -30*dB
        instrument.mwb.config[qubit_name].c1atten         = 31*dB
        instrument.mwb.config[qubit_name].c2gain          = -30*dB
        instrument.mwb.config[qubit_name].c2atten         = 31*dB
        instrument.qla.config[qubit_name].window          = (note.cavity_readout_window_weight_I, note.cavity_readout_window_weight_Q)
        instrument.qla.status.cooltime                    = 50*us
        instrument.qla.set_acquisition_mode(averaging_shots = False, averaging_waveform = True)
        instrument.seq[qubit_name].readout.delay          = note.cavity_readout_trigger_delay + note.cavity_readout_window_delay
        instrument.seq[qubit_name].readout.duration       = note.cavity_readout_window_length

        instrument.seq[qubit_name].qubit.seq              = "T"
        instrument.seq[qubit_name].readout.seq            = "T"
        instrument.seq[qubit_name].cr1.seq                = "T"
        instrument.seq[qubit_name].cr2.seq                = "T"

        if qubit_name in qubit_notes.keys():
            alpha = name_to_alpha(qubit_name)
            instrument.seq.set_user_command("HPI{0}".format(alpha), "A{0} DT{1} B{2} D{2} B{2}".format(note.pi_pulse_power, note.half_pi_pulse_drag_coeff, note.half_pi_pulse_length_precise["ns"]))
            instrument.seq.set_user_command("MEAS{0}".format(alpha), "A{0} T B10 M F3000".format(note.cavity_readout_window_power))
            instrument.seq.set_user_command("WAIT{0}".format(alpha), "B{0} ".format(note.half_pi_pulse_length_precise["ns"]))

def set_cross(instrument, cross_notes, cross_information):
    """Setup of the instruments used in the experiment
    Args:
        instrument (TimeDomainInstrumentManager): Instruments used in the experiment
        cross_notes (AttributionDict): Calibration notes for two-qubit control
        cross_information (list): List of the cross-resonance port names
    """
    for cross_name in cross_information:
        if str(cross_name) in cross_notes.keys():
            note = cross_notes[str(cross_name)]
            if cross_name[2] is "cr1":
                instrument.mwb.config[cross_name[1]].c1gain          = note.gain
                instrument.mwb.config[cross_name[1]].c1atten         = note.atten
            if cross_name[2] is "cr2":
                instrument.mwb.config[cross_name[1]].c2gain          = note.gain
                instrument.mwb.config[cross_name[1]].c2atten         = note.atten
            calpha = name_to_alpha(cross_name[0])
            talpha = name_to_alpha(cross_name[1])
            instrument.seq.set_user_command("DCR{0}{1}".format(calpha, talpha), "B{0} FT{1} A{2} P{3} F{4} B{0}".format(note.crw["ns"], note.crr["ns"], note.cra, note.crp, note.crt["ns"]))
            instrument.seq.set_user_command("DCT{0}{1}".format(calpha, talpha), "B{0} FT{1} A{2} P{3} F{4} B{0}".format(note.crw["ns"], note.crr["ns"], note.cta, note.ctp, note.crt["ns"]))
            instrument.seq.set_user_command("WCR{0}{1}".format(calpha, talpha), "B{0}                 B{1} B{0}".format(note.crw["ns"], note.crt["ns"]))
            instrument.seq.set_user_command("WAIT{0}{1}".format(calpha, talpha), "B{0} ".format(note.crw["ns"] + 0.5*note.crt["ns"]))

@retry(stop_max_attempt_number=5, wait_fixed=60)
def take_data(job_table):
    instrument = tdm_inst
    target_list = qubit_information
    exp       = CreateProjectorTarget(target_list)
    projector = exp.execute(tdm_inst, qubit_notes, save=False)
    set_qubit(instrument, qubit_notes, qubit_information)
    set_cross(instrument, cross_notes, cross_information)
    port_list = job_table.table[0].sequence.keys()
    
    sequences = {}
    for port_name in port_list:
        if type(port_name) is str:
            sequences[(port_name, "readout")] = []
            sequences[(port_name,   "qubit")] = []
            sequences[(port_name,     "cr1")] = []
            sequences[(port_name,     "cr2")] = []
    for job in job_table.table:
        for port_name in job.sequence.keys():
            if type(port_name) is str:
                alpha = name_to_alpha(port_name)
                sequences[(port_name, "readout")].append("MEAS{0}".format(alpha))
                sequences[(port_name,   "qubit")].append(job.sequence[port_name] + " T")
                sequences[(port_name,     "cr1")].append("T")
                sequences[(port_name,     "cr2")].append("T")
            if type(port_name) is tuple:
                if port_name[1] in target_list:
                    del sequences[(port_name[1], port_name[2])][-1]
                    sequences[(port_name[1], port_name[2])].append(job.sequence[port_name] + " T")
                    
    sweep_index = 0
    sweep_axis = []
    for port_name in port_list:
        if type(port_name) is str:
            instrument.seq.config_variable_command("C{0}C".format(sweep_index+0), sequences[(port_name, "readout")], "readout_sequence")
            instrument.seq.config_variable_command("C{0}C".format(sweep_index+1), sequences[(port_name,   "qubit")],   "qubit_sequence")
            instrument.seq.config_variable_command("C{0}C".format(sweep_index+2), sequences[(port_name,     "cr1")],     "cr1_sequence")
            instrument.seq.config_variable_command("C{0}C".format(sweep_index+3), sequences[(port_name,     "cr2")],     "cr2_sequence")
            instrument.seq[port_name].readout.seq = "C{0}C".format(sweep_index+0)
            instrument.seq[port_name].qubit.seq   = "C{0}C".format(sweep_index+1)
            instrument.seq[port_name].cr1.seq     = "C{0}C".format(sweep_index+2)
            instrument.seq[port_name].cr2.seq     = "C{0}C".format(sweep_index+3)
            sweep_index += 4
            sweep_axis += [0]*4
            
    instrument.qla.status.shots = 1024
    dataset = instrument.take_data(dataset_name="test",save=False,sweep_axis=sweep_axis)
    data_dict = dataset.get_iq_data_dict()
    hist_dict = projector.get_histogram_dict(data_dict)
    for job, histogram in zip(job_table.table, hist_dict):
        job.result = normalize_histogram(histogram)