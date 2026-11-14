"""Independent face-frame physical rates from a public-input snapshot."""
import copy,math
import numpy as np
import optical_check as oc

def _strengths(catalog, operator, stability_tensor, frames):
    """Excitation quadratic form in the specimen-face frame, not lab frame."""
    result = np.zeros((len(catalog), 72), dtype=float)
    truth = operator.TRUTH
    instrument = operator.INSTRUMENT
    base = copy.deepcopy(truth["surface_stack"])
    try:
        for k, (target, _, laser, azimuth, tx, ty, plate, _) in enumerate(catalog):
            if int(target) == 0:
                continue
            laser = int(laser)
            n = instrument["sample_refractive_index"][str(laser)]
            mirror = instrument["dichroic_reflection_at_laser_line"][str(laser)]
            polarization = math.radians(2.0 * plate)
            phase = math.radians(mirror["phase_s_minus_p_deg"])
            incoming = np.array([
                math.sqrt(mirror["R_p"]) * math.cos(polarization),
                math.sqrt(mirror["R_s"]) * complex(math.cos(phase), math.sin(phase))
                * math.sin(polarization), 0.0])
            for branch in truth["optical_states"]["branches"]:
                b = int(branch["branch_id"])
                thickness = branch["common_layer_thickness_scale"]
                truth["surface_stack"]["sample_layers"] = [
                    dict(layer, thickness_nm=thickness * layer["thickness_nm"])
                    for layer in base["sample_layers"]]
                stage = (operator._turn("y", ty + branch["tilt_y_registration_deg"]) @ operator._turn("x", tx + branch["tilt_x_registration_deg"])
                         @ operator._turn("z", azimuth + branch["azimuth_registration_deg"]))
                field = operator._field_in_face_frame(incoming, stage, n, laser, "sample")
                for d, frame in enumerate(frames):
                    in_face = frame @ stability_tensor @ frame.T
                    result[k, 6*d+b] = float(np.real(
                        field.conjugate() @ in_face @ field) / (2.0 / (1.0+n))**2)
    finally:
        truth["surface_stack"] = base
    return result


def independent_rates():
    model=oc.TRUTH['photodynamics'];catalog,intact,instrument=oc.independent_rates();cfg=copy.deepcopy(instrument)
    A=np.asarray(model['local_strength_tensor']);j0=_strengths(catalog,oc,A,oc.VARIANTS)
    unit=cfg['exposure_unit_s'];c=model['specimen_collection_scale'];s=model['specimen_signal_scale']
    result=np.repeat(np.repeat(intact[:,:,None,:]/unit,3,axis=2),2,axis=1);Q=np.zeros((len(catalog),144,3,3))
    original=oc.VARIANTS
    try:
        for k,family in enumerate(model['kinetic_families']):
            rotation=np.array(family['switched_parent_rotation'],float)
            altered=[rotation@frame for frame in original];oc.VARIANTS=altered
            _,switched,_=oc.independent_rates();j1=_strengths(catalog,oc,A,altered)
            for i,row in enumerate(catalog):
                if int(row[0])==0:continue
                laser=str(int(row[2]));bg=np.array(cfg['sample_background_electrons_per_exposure_unit'][laser])
                bright=c*(bg+s*(intact[i]-bg));changed=c*(bg+model['dark_signal_fraction'][laser]*s*(switched[i]-bg))
                average=np.tile(bright.reshape(12,6,4).mean(axis=0),(12,1))
                result[i,k::2,0]=bright/unit;result[i,k::2,1]=changed/unit;result[i,k::2,2]=average/unit
                co=family['transition_coefficients_per_s_per_strength'][laser]
                for a,b,key,field in ((0,1,'bright_to_dark',j0),(0,2,'bright_to_scrambled',j0),(1,0,'dark_to_bright',j1),(1,2,'dark_to_scrambled',j1)):
                    Q[i,k::2,a,b]=field[i]*co[key]
                for a in (0,1):Q[i,k::2,a,a]=-Q[i,k::2,a].sum(axis=-1)
    finally:oc.VARIANTS=original
    return catalog,result,Q,cfg

def build():return independent_rates()
