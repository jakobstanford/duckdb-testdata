create TABLE  calcs AS (
select
db.*,


/* === gather demands (Qud and Quf) === */
quf.object_t_value, quf.element_t_value,
quf.step_type AS step_type,
quf.case_combo,
quf.result_case_name AS result_case_name_quf,
qud.result_case_name AS result_case_name_qud,


/* = considerations for modelling of tension-only rods =*/
sign(qud.axial) * least(abs(qud.axial) + P_rod_delta, abs(qud.axial) * P_rod_factor) AS P_UD,
sign(quf.axial) * least(abs(quf.axial) + P_rod_delta, abs(quf.axial) * P_rod_factor) AS P_UF,

/* === determine demands === */

/* additional bending moment allowance for load-path connection eccentricity  (ref. detail 7/S8.13) */
CASE
 WHEN db.beam_name = '1225' THEN qud.mx + 6000*sign(qud.mx)
 WHEN db.beam_name IN ('791', '793') THEN qud.mx + 7500*sign(qud.mx)
 WHEN db.beam_name IN ('789', '775', '797', '795') THEN qud.mx + 4800*sign(qud.mx)
 WHEN db.beam_name = '481' and is_seismic_combo THEN qud.mx + 10224*sign(qud.mx)
 WHEN db.beam_name = '481' and not is_seismic_combo THEN qud.mx + 5400*sign(qud.mx)
 ELSE qud.mx END AS M_UDx,
CASE
 WHEN db.beam_name = '1225' THEN quf.mx + 6000*sign(quf.mx)
 WHEN db.beam_name IN ('791', '793') THEN quf.mx + 7500*sign(quf.mx)
 WHEN db.beam_name IN ('789', '775', '797', '795') THEN quf.mx + 4800*sign(quf.mx)
 WHEN db.beam_name = '481' and is_seismic_combo THEN quf.mx + 10224*sign(quf.mx)
 WHEN db.beam_name = '481' and not is_seismic_combo THEN quf.mx + 5400*sign(quf.mx)
 ELSE quf.mx END AS M_UFx,

qud.my AS M_UDy,
quf.my AS M_UFy,

qud.torque AS M_UDt,
quf.torque AS M_UFt,

qud.vx AS V_UDx,
quf.vx AS V_UFx,

/* additional 300k allowance for load-path connection eccentricity  (ref. detail 7/S8.13) */
CASE WHEN db.beam_name = '1225' THEN qud.vy + 300*sign(qud.vy) ELSE qud.vy END AS V_UDy,
CASE WHEN db.beam_name = '1225' THEN quf.vy + 300*sign(quf.vy) ELSE quf.vy END AS V_UFy,



/* === determine capacities (Qce and Qcl) === */


/* phi correction factors for Q_cl:
 - assume Q_cl was computed with phi per AISC360.
 - for ASCE41, phi=1 though, hence divide capacities by phi.
 - this applies only in case of seismic (ASCE41) combos. for ASCE7 combos, do not apply this correction (set correction factor to 1) */

True AS allow_phi_1,  /* toggle whether phi _may_ be =1 for force controlled cases under seismic-load combos (non-seismic combos are treated force-controlled and use capacities with phi<1. */

CASE WHEN db.prop_name LIKE 'W%' AND db.prop_name NOT LIKE 'WT%' THEN
  (ep_d - 2*ep_tf - 2*ep_rf) / ep_tw
ELSE
  9999.0   /* large value to not trigger Eq G2-2*/
END AS h_tw,  /* h/t_w as used in AISC Eq G2-2*/

CASE WHEN is_seismic_combo AND allow_phi_1 THEN 0.9 ELSE 1 END AS phi_qcl_corr,  /* phi correction for bending&axial  */

CASE WHEN is_seismic_combo AND allow_phi_1
  THEN  /* check short web (AISC360, G2.1(a)) */
    CASE WHEN db.prop_name LIKE 'W%' AND db.prop_name NOT LIKE 'WT%' AND h_tw < 2.24 * sqrt(29000 / fy) THEN 1.0 ELSE 0.9 END  /* apply either correction factor */
  ELSE
   1 /* (no correction factor) */
END AS phi_qcl_corr_vy,  /* phi correction for shear in Y direction  */

CASE WHEN is_seismic_combo AND allow_phi_1 THEN 0.9 ELSE 1.0 END AS phi_qcl_corr_vx,  /* phi correction for shear in X direction  */

CASE WHEN is_seismic_combo AND allow_phi_1 THEN 1 ELSE 0.75 END AS phi_qcl_tens_rupture,  /* phi factor for tensile rupture (not a correction factor) */


qcl_pccomp / phi_qcl_corr AS Pcc_cl, /* lower bound compressive strength */


CASE WHEN member_type_asce41 = 'brace'
 THEN least(qcl_pctension, phi_qcl_tens_rupture * An*1*58)
 ELSE qcl_pctension / phi_qcl_corr
END AS Pct_CL, /* lower bound tensile strength */

CASE WHEN Pct_CL < qcl_pctension  THEN 'net_rupt'  ELSE 'gross_yld' END AS Pct_CL_gov_by,

qce_pccomp AS Pcc_CE, /* expected compressive strength */

CASE WHEN member_type_asce41 = 'brace'
 THEN least(qce_pctension, 1.0*An*1*(1.2*58))  /* assumed U=1, phi=1 */
 ELSE qce_pctension
END AS Pct_CE, /* expected tensile strength */

CASE WHEN Pct_CE < qce_pctension THEN 'net_rupt' ELSE 'gross_yld' END AS Pct_CE_gov_by,

CASE WHEN (member_type_asce41 = 'brace' AND is_seismic_combo) THEN Pcc_CE ELSE Pcc_CL END AS Pcc,
CASE WHEN (member_type_asce41 = 'brace' AND is_seismic_combo) THEN Pct_CE ELSE Pct_CL END AS Pct,

CASE WHEN P_UF < 0 THEN Pcc ELSE Pct END AS Pc_quf,
CASE WHEN P_UD < 0 THEN Pcc ELSE Pct END AS Pc_qud,


qce_mcmajor AS Mcex_AISC,
qce_mcminor AS Mcey,
qcl_mcmajor / phi_qcl_corr AS Mclx,
qcl_mcminor / phi_qcl_corr AS Mcly,

/*check if flexural strength Mcex is governed by shear yielding of web*/
CASE
  WHEN member_type_asce41 = 'column' THEN 0.6 * Fye * (ep_d - (2 * ep_tf))*ep_tw
  ELSE NULL
END AS Vce_web,

CASE
  WHEN (Vce_web IS NOT NULL AND member_flex_releases = 'Fix-Fix') THEN Vce_web * (L/2.0)
  WHEN (Vce_web IS NOT NULL AND member_flex_releases = 'Fix-Pin') THEN Vce_web * L
  WHEN (Vce_web IS NOT NULL AND member_flex_releases = 'Pin-Fix') THEN Vce_web * L
  WHEN (Vce_web IS NOT NULL AND member_flex_releases = 'Pin-Pin') THEN qce_mcmajor*1000
  WHEN (Vce_web IS NOT NULL AND member_flex_releases = 'Other'  ) THEN 0
  ELSE NULL
END AS Mce_web,

Mce_web / Mcex_AISC AS Mcew_Mcex,

CASE
  WHEN Mcew_Mcex < 1 THEN Mce_web
  ELSE Mcex_AISC
END AS Mcex,

CASE
  WHEN Mcew_Mcex < 1 THEN 'Shear Yielding of Web'
  ELSE 'AISC Chapter F'
END AS Mcex_gov_by,

qces_phivnminor AS Vcex,
qces_phivnmajor AS Vcey,
qcls_phivnminor / phi_qcl_corr_vx AS Vclx,
qcls_phivnmajor / phi_qcl_corr_vy AS Vcly,

abs(P_UF) / Pye AS Puf_Pye,

/* === determine initial m-factors === */
CASE

  WHEN member_type LIKE '%-diag' THEN 2.0

  WHEN member_type LIKE '%-vert' AND db.prop_name NOT IN ('W14X90','W14X109') THEN 1.5
  WHEN member_type LIKE '%-vert' AND db.prop_name IN     ('W14X90','W14X109') THEN
    CASE WHEN Puf_Pye < 0.2 THEN 1.25 ELSE 1.0  END

  WHEN member_type LIKE '%-horiz' AND db.prop_name NOT IN ('W14X26','W14X30','W14X34','W14X43','W14X90','W14X109','W18X35','W18X76') THEN 2.0
  WHEN member_type LIKE '%-horiz' AND db.prop_name IN     ('W14X26','W14X30','W14X34','W14X43','W14X90','W14X109','W18X35','W18X76') THEN
    CASE WHEN Puf_Pye < 0.2 THEN 1.25 ELSE 1.0  END

  WHEN member_type LIKE '%-rod' THEN 2.0
END AS m_initial,



/* === check plastic hinge (9-10, 9-11) === */
/* do not check plastic hinge if Puf_Pye > 0.6 */
/* OK to use MUDy and m_y_hinge > 1 when Mcex<Mpex and is_seismic_combo*/
ep_Zx * Fye AS Mpex,
ep_Zy * Fye AS Mpey,
ep_Zx * fy  AS Mplx,

Mcex / Mpex AS mcex_mpex,  /* Mcex / Mpex < 1 is same as Mcex < Mpex */

CASE
  WHEN member_type_asce41 = 'brace' THEN NULL
  WHEN Puf_Pye > 0.6 OR NOT is_seismic_combo OR m_initial = 1.0 THEN NULL
  WHEN mcex_mpex < 0.999 AND Puf_Pye <= 0.6 THEN 'qufx_qudy'
  ELSE 'qud'
END AS hinge_check_forces,


CASE
  WHEN hinge_check_forces = 'qud'    THEN m_initial
  WHEN hinge_check_forces = 'qufx_qudy' THEN 1.0
  ELSE NULL
END AS m_x_hinge,

CASE
  WHEN hinge_check_forces IS NOT NULL THEN m_initial
  ELSE NULL
END AS m_y_hinge,


CASE
  WHEN hinge_check_forces = 'qud' THEN M_UDx
  WHEN hinge_check_forces = 'qufx_qudy' THEN M_UFx
  ELSE NULL
END AS Mux_hinge,

CASE
  WHEN hinge_check_forces IS NOT NULL THEN M_UDy
  ELSE NULL
END AS Muy_hinge,


CASE
  WHEN hinge_check_forces = 'qud' THEN Mpex
  WHEN hinge_check_forces = 'qufx_qudy' THEN Mplx
  ELSE NULL
END AS Mpx_hinge,

CASE
  WHEN hinge_check_forces IS NOT NULL THEN Mpey
  ELSE NULL
END AS Mpy_hinge,


CASE
  WHEN hinge_check_forces IS NOT NULL THEN
    CASE WHEN Puf_Pye < 0.2 * kappa THEN
     (  Puf_Pye / 2.0 +             ( abs(Mux_hinge) / (m_x_hinge * Mpx_hinge)  +  abs(Muy_hinge) / (m_y_hinge * Mpy_hinge)) ) / kappa
    ELSE
     (  Puf_Pye       + .88889::float * ( abs(Mux_hinge) / (m_x_hinge * Mpx_hinge)  +  abs(Muy_hinge) / (m_y_hinge * Mpy_hinge)) ) / kappa
    END

  ELSE NULL
END AS dcr_hinge,



/* === check 'column' member (9-12, 9-13, 9-14, 9-15) === */
/* = PMM interaction =*/
CASE WHEN (Puf_Pye > 0.6 AND member_type_asce41 = 'column') OR NOT is_seismic_combo OR m_initial = 1.0 THEN 'quf' ELSE 'qud' END AS mem_check_forces,

CASE WHEN mem_check_forces = 'quf' THEN 1.0 ELSE m_initial END AS m_mem,
CASE WHEN member_type_asce41 = 'column' THEN 1.0 ELSE m_initial END AS m_P_mem_flex,  /*TODO: add below*/

CASE WHEN mem_check_forces = 'quf' THEN M_UFx ELSE M_UDx END AS Mux_mem,
CASE WHEN mem_check_forces = 'quf' THEN M_UFy ELSE M_UDy END AS Muy_mem,

CASE
  WHEN (member_type_asce41 = 'brace' AND mem_check_forces = 'qud') THEN P_UD
  ELSE P_UF
END AS Pu_mem_flex,

CASE WHEN mem_check_forces = 'quf' THEN Mclx ELSE Mcex END AS Mcx_mem,
CASE WHEN mem_check_forces = 'quf' THEN Mcly ELSE Mcey END AS Mcy_mem,


CASE
  WHEN (member_type_asce41 = 'brace' AND mem_check_forces = 'qud') THEN Pc_qud
  ELSE Pc_quf
END AS Pc_mem_flex,


abs(Pu_mem_flex) / (m_P_mem_flex * Pc_mem_flex) AS dcr_mem_flex_p,
abs(Mux_mem) / (m_mem * Mcx_mem) AS dcr_mem_flex_mx,
abs(Muy_mem) / (m_mem * Mcy_mem) AS dcr_mem_flex_my,
CASE WHEN rt_dc IS NOT NULL THEN rt_dc ELSE 0 END AS dcr_mem_flex_mt,

CASE WHEN member_type NOT LIKE '%-rod' THEN
  CASE WHEN dcr_mem_flex_p < 0.2 * kappa THEN
   ( (dcr_mem_flex_p / 2.) +             ( dcr_mem_flex_mx  +  dcr_mem_flex_my + dcr_mem_flex_mt ) ) / kappa
  ELSE
   ( dcr_mem_flex_p       + .88889::float * ( dcr_mem_flex_mx  +  dcr_mem_flex_my + dcr_mem_flex_mt ) ) / kappa
  END
END AS dcr_mem_flex,


/* = axial check =*/
CASE
  WHEN member_type_asce41 = 'column' AND is_seismic_combo THEN
    CASE WHEN P_UF > 0 THEN
      (abs(P_UD) / (m_mem * Pc_qud)) / kappa   /* tension - (9-15) */
    ELSE
      Puf_Pye / (0.75 * kappa)   /* compression - (9-14) */
    END
  WHEN member_type_asce41 IN ('brace','column') AND NOT is_seismic_combo AND P_UF < 0 THEN
    abs(P_UF) / (Pcc_CL * 0.75 * kappa)
  ELSE
    NULL /* otherwise check is not applicable */
END AS dcr_mem_axial,

/* = rod checks = */
CASE
  WHEN db.prop_name = '1inROD' THEN 43.2
  WHEN db.prop_name = '1.25inROD' THEN 70.0
  WHEN db.prop_name = '1.625inROD' THEN 135.9
END AS Tye_rod,

CASE WHEN member_type LIKE '%-rod' THEN 2 * abs(P_ud) / (m_initial * Tye_rod) END AS dcr_mem_rod,


greatest(dcr_mem_axial, dcr_mem_flex, dcr_mem_rod) AS dcr_mem,



/* === shear checks === */
CASE WHEN mem_check_forces = 'quf' THEN V_UFx ELSE V_UDx END AS Vux_mem,
CASE WHEN mem_check_forces = 'quf' THEN V_UFy ELSE V_UDy END AS Vuy_mem,
CASE WHEN mem_check_forces = 'quf' THEN Vclx ELSE Vcex END AS Vcx,
CASE WHEN mem_check_forces = 'quf' THEN Vcly ELSE Vcey END AS Vcy,

abs(Vux_mem) / Vcx AS dcr_shear_x,
abs(Vuy_mem) / Vcy AS dcr_shear_y,
greatest(dcr_shear_x, dcr_shear_y) AS dcr_shear,



/* === DCR on the main "deformation controlled" criterion === */
CASE WHEN member_type_asce41 = 'column' THEN
  Puf_Pye / 0.60
ELSE
  NULL
END AS dcr_defocntl_check,


/* === determine governing DCR === */
greatest(dcr_hinge, dcr_mem,  dcr_defocntl_check, dcr_shear_x, dcr_shear_y ) AS dcr,

CASE
 WHEN dcr = dcr_hinge THEN 'hinge'
 WHEN dcr = dcr_defocntl_check THEN 'defocntl'
 WHEN dcr = dcr_mem THEN 'member'
 WHEN dcr = dcr_shear_x THEN 'shear x'
 WHEN dcr = dcr_shear_y THEN 'shear y'
END AS dcr_gvt_by,


/* === row numbering to filter out governing results === */
row_number() over (partition BY quf.beam_name ORDER BY dcr DESC) AS rn

FROM beam_list AS db

JOIN results_quf AS quf ON (db.beam_name = quf.beam_name)

LEFT JOIN results_qud AS qud ON (
  quf.beam_name = qud.beam_name AND
  quf.object_t_value = qud.object_t_value AND
  quf.element_t_value = qud.element_t_value AND
  quf.case_combo = qud.case_combo AND
  quf.step_type = qud.step_type
  )
);
