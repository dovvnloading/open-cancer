"""Invented-data checks of prediction isolation and intervention semantics."""
from __future__ import annotations

import hashlib
import importlib.util
import itertools
import math
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
import numpy as np


HERE = Path(__file__).resolve().parents[1]/"research/barcode-lineage-memory"


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE/filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


MODEL = load("barcode_model_test", "model.py")
with patch.dict(sys.modules, {"model": MODEL}):
    ANALYSIS = load("barcode_analysis_test", "analysis.py")
with patch.dict(sys.modules, {"model": MODEL, "analysis": ANALYSIS}):
    EXPLORATORY = load("barcode_exploratory_test", "exploratory.py")


def invented_counts():
    p1 = np.array([[60,20,10,5,5], [10,60,20,5,5], [20,10,60,5,5], [30,30,30,5,5]])
    p2 = np.array([[12,15,63,5,5], [50,10,30,5,5], [45,45,0,5,5], [1,2,87,5,5]])
    counts = np.empty((5,8), dtype=np.int64)
    counts[:, ::2], counts[:, 1::2] = p1.T, p2.T
    return counts


def reference_nested(early, late, target, delta=1e-4, donor_early=None):
    """Small probability-space enumeration with standard-library summation."""
    def unit(values):
        total = math.fsum(values)
        return [float(value)/total for value in values]
    def smoothed(values):
        p = unit(values)
        return [(1-delta)*value+delta/len(p) for value in p]
    x = [smoothed(row) for row in early]
    d = x if donor_early is None else [smoothed(row) for row in donor_early]
    y = [unit(row) for row in late]
    training = [r for r in range(4) if r != target]
    def forecast(row, donors, a, b):
        z = [math.fsum(d[r][i] for r in donors)/len(donors) for i in range(len(x[row]))]
        return unit([own**a*other**b for own,other in zip(x[row], z)])
    candidates = []
    for a,b in itertools.product((.5,1,1.5,2,3,4), (0,.25,.5,1,2)):
        scores=[]
        for inner in training:
            q=forecast(inner, [r for r in training if r != inner], a,b)
            scores.append(math.fsum(p*math.log(pred) for p,pred in zip(y[inner],q)))
        candidates.append((a,b,math.fsum(scores)/3))
    def select(rows):
        best=max(row[2] for row in rows)
        return min((row for row in rows if best-row[2]<=1e-12),key=lambda row:(row[1],abs(row[0]-1),row[0]))
    a0,_,inner0=select([row for row in candidates if row[1]==0])
    a,b,inner1=select(candidates)
    q0,q1=forecast(target,training,a0,0),forecast(target,training,a,b)
    gain=math.fsum(p*math.log(full/base) for p,full,base in zip(y[target],q1,q0))
    return a0,a,b,inner0,inner1,np.asarray(q0),np.asarray(q1),gain


def temporal_fixture():
    # Forty artificial identifiers split evenly by a hash criterion. Frequencies
    # create one home and three minor receiver profiles, with nonempty halves.
    values=[]
    desired=[0,0]
    for i in range(1000):
        name=f"SYNTHETIC_BARCODE_{i:06d}"
        half=hashlib.sha256(name.encode()).digest()[-1]%2
        if desired[half]<20:
            values.append(name);desired[half]+=1
        if desired==[20,20]:
            break
    values.append("SYNTHETIC_BACKGROUND")
    raw=np.ones((41,8),dtype=np.int64)
    raw[:40,0]=np.arange(1,41)*100
    raw[:40,1]=np.arange(40,0,-1)*100
    raw[40,0:2]=100
    for r in range(1,4):
        raw[:40,2*r]=np.arange(40,0,-1)
        raw[40,2*r]=10**7
        raw[40,2*r+1]=10**7
    return raw,np.asarray(values)


class BarcodeForecastTests(unittest.TestCase):
    def test_all_folds_match_independent_probability_enumeration(self):
        obs=MODEL.prepare_observations(invented_counts())
        for delta,target in itertools.product((1e-6,1e-4,1e-3),range(4)):
            fit=MODEL.fit_outer_fold(obs.early,obs.late,target,delta)
            a0,a,b,s0,s1,q0,q1,gain=reference_nested(obs.early,obs.late,target,delta)
            self.assertEqual((fit.base_a,fit.full_a,fit.full_b),(a0,a,b))
            self.assertAlmostEqual(fit.base_inner_score,s0,places=12)
            self.assertAlmostEqual(fit.full_inner_score,s1,places=12)
            np.testing.assert_allclose(np.exp(fit.base_log_prediction),q0,rtol=1e-12)
            np.testing.assert_allclose(np.exp(fit.full_log_prediction),q1,rtol=1e-12)
            self.assertAlmostEqual(fit.gain,gain,places=12)

    def test_held_out_late_cannot_change_prediction(self):
        obs=MODEL.prepare_observations(invented_counts())
        for target in range(4):
            altered=obs.late.copy()
            altered[target]=[1,0,0,0,0,0]
            before=MODEL.fit_outer_fold(obs.early,obs.late,target)
            after=MODEL.fit_outer_fold(obs.early,altered,target)
            self.assertEqual((before.base_a,before.full_a,before.full_b),(after.base_a,after.full_a,after.full_b))
            np.testing.assert_array_equal(before.full_log_prediction,after.full_log_prediction)
            self.assertNotEqual(before.gain,after.gain)

    def test_outer_early_excluded_from_inner_parameter_choice(self):
        obs=MODEL.prepare_observations(invented_counts())
        for target in range(4):
            changed=obs.early.copy();changed[target]=[0,0,0,0,1,0]
            before=MODEL.fit_outer_fold(obs.early,obs.late,target)
            after=MODEL.fit_outer_fold(changed,obs.late,target)
            self.assertEqual((before.base_a,before.full_a,before.full_b),(after.base_a,after.full_a,after.full_b))
            self.assertEqual(before.base_inner_score,after.base_inner_score)
            self.assertEqual(before.full_inner_score,after.full_inner_score)

    def test_p1_only_dictionary_and_other_unseen(self):
        raw=np.vstack([invented_counts(),[0,10,0,20,0,30,0,40]])
        obs=MODEL.prepare_observations(raw)
        self.assertEqual(obs.early_dictionary_lineages,5)
        np.testing.assert_array_equal(obs.early[:,-1],np.zeros(4))
        np.testing.assert_allclose(obs.late[:,-1],[10/110,20/120,30/130,40/140])
        raw[-1,1::2]*=100
        altered=MODEL.prepare_observations(raw)
        np.testing.assert_array_equal(obs.keep_mask,altered.keep_mask)
        np.testing.assert_array_equal(obs.early,altered.early)

    def test_threshold_keeps_boundary_and_recomputes_support(self):
        raw=np.tile(np.array([50,25,15,10]),(8,1)).T
        obs=MODEL.prepare_observations(raw,.15)
        self.assertEqual(obs.early_dictionary_lineages,3)
        np.testing.assert_allclose(obs.retained_mass,.9)
        np.testing.assert_allclose(obs.early[0],[50/90,25/90,15/90,0])

    def test_empty_samples_and_invalid_inputs_are_not_dropped(self):
        with self.assertRaises(MODEL.EmptySampleError): MODEL.prepare_observations(invented_counts(),.99)
        for value in (-.1,1,float("nan")):
            with self.assertRaises(ValueError): MODEL.prepare_observations(invented_counts(),value)
        for value in (-1,.1,float("nan")):
            raw=invented_counts().astype(float);raw[0,0]=value
            with self.assertRaises(ValueError): MODEL.prepare_observations(raw)
        with self.assertRaises(ValueError): MODEL.prepare_observations(np.zeros((1,7)))

    def test_tie_breaking_prioritizes_lower_b_then_a_near_one(self):
        rows=[(4,0,0),(.5,0,0),(1.5,0,0),(1,.25,0),(1,0,-5e-13)]
        self.assertEqual(MODEL.select_candidate(rows),(1,0,-5e-13))
        self.assertEqual(MODEL.select_candidate([(.5,0,0),(1.5,0,0)]),(.5,0,0))

    def test_culture_read_depth_does_not_change_nested_fit(self):
        raw=invented_counts()
        obs=MODEL.prepare_observations(raw)
        altered=MODEL.prepare_observations(raw*np.array([1,1000,2,100,3,10,4,1]))
        for target in range(4):
            first=MODEL.fit_outer_fold(obs.early,obs.late,target)
            second=MODEL.fit_outer_fold(altered.early,altered.late,target)
            self.assertAlmostEqual(first.gain,second.gain,places=13)
            self.assertEqual((first.base_a,first.full_a,first.full_b),(second.base_a,second.full_a,second.full_b))

    def test_donor_only_intervention_preserves_baseline_and_observations(self):
        raw=invented_counts()
        original=MODEL.prepare_observations(raw)
        data=EXPLORATORY.fixed_dictionary_intervention(raw,.06,"donor_feature_only")
        np.testing.assert_array_equal(data["early"],original.early)
        np.testing.assert_array_equal(data["late"],original.late)
        self.assertEqual(data["early_dictionary_lineages"],original.early_dictionary_lineages)
        for target in range(4):
            baseline=MODEL.fit_outer_fold(original.early,original.late,target)
            new=MODEL.fit_outer_fold(data["early"],data["late"],target,donor_early=data["donor_early"])
            self.assertEqual(baseline.base_a,new.base_a)
            np.testing.assert_array_equal(baseline.base_log_prediction,new.base_log_prediction)
            expected=reference_nested(data["early"],data["late"],target,donor_early=data["donor_early"])
            self.assertAlmostEqual(new.gain,expected[-1],places=12)

    def test_late_threshold_precedes_other_category_collapse(self):
        # The two late-only lineages individually fall below .02, although their
        # collapsed mass exceeds it. Both must be removed before OTHER collapse.
        raw=np.tile(np.array([100,0,0]),(8,1)).T
        raw[1:,1::2]=2
        data=EXPLORATORY.fixed_dictionary_intervention(raw,.02,"late_only")
        self.assertEqual(data["early"].shape,(4,2))
        np.testing.assert_allclose(data["late"],[[1,0]]*4)
        np.testing.assert_allclose(data["retained_mass"][:,1],100/104)

    def test_fixed_dictionary_survives_early_threshold(self):
        data=EXPLORATORY.fixed_dictionary_intervention(invented_counts(),.06,"early_only")
        self.assertEqual(data["early"].shape,(4,6))
        np.testing.assert_array_equal(data["early"][:,3:],np.zeros((4,3)))
        self.assertEqual(data["early_dictionary_lineages"],5)

    def test_probability_smoothing_and_forecasts_normalize(self):
        obs=MODEL.prepare_observations(invented_counts())
        fit=MODEL.fit_outer_fold(obs.early,obs.late,0)
        for logs in (fit.base_log_prediction,fit.full_log_prediction):
            self.assertTrue(np.isfinite(logs).all())
            self.assertAlmostEqual(float(np.exp(logs).sum()),1,places=13)
        for delta in (0,1,-.01,float("nan")):
            with self.assertRaises(ValueError): MODEL.fit_outer_fold(obs.early,obs.late,0,delta)

    def test_temporal_hash_half_assignment_is_exact(self):
        raw,barcodes=temporal_fixture()
        expected=[hashlib.sha256(value.encode("ascii")).digest()[-1]&1 for value in barcodes]
        self.assertEqual(EXPLORATORY.half_parity(barcodes).tolist(),expected)

    def test_temporal_selected_support_does_not_depend_on_late_counts(self):
        raw,barcodes=temporal_fixture()
        first=EXPLORATORY.temporal_diagnostic(raw,barcodes)
        altered=raw.copy();altered[:,1::2]*=3
        second=EXPLORATORY.temporal_diagnostic(altered,barcodes)
        for left,right in zip(first["home_results"],second["home_results"]):
            self.assertEqual(left["selected_barcodes"],right["selected_barcodes"])
            self.assertEqual(left["half_selected_barcodes"],right["half_selected_barcodes"])

    def test_temporal_test_receiver_half_never_selects_its_weight(self):
        raw,barcodes=temporal_fixture()
        first=EXPLORATORY.temporal_diagnostic(raw,barcodes)["home_results"][0]
        self.assertTrue(first["available"])
        parity=EXPLORATORY.half_parity(barcodes[:40])
        changed=raw.copy();changed[np.flatnonzero(parity==0),2]=1
        second=EXPLORATORY.temporal_diagnostic(changed,barcodes)["home_results"][0]
        one=next(r for r in first["receiver_half_results"] if r["receiver_culture"]=="DT2" and r["test_half"]==0)
        two=next(r for r in second["receiver_half_results"] if r["receiver_culture"]=="DT2" and r["test_half"]==0)
        self.assertEqual(one["selected_weight"],two["selected_weight"])
        self.assertEqual(one["receiver_train_read_count"],two["receiver_train_read_count"])

    def test_temporal_sparse_support_is_counted_not_invented(self):
        raw,barcodes=temporal_fixture()
        parity=EXPLORATORY.half_parity(barcodes[:40])
        raw[:40,2]=0
        for half in (0,1): raw[np.flatnonzero(parity==half)[0],2]=5
        home=EXPLORATORY.temporal_diagnostic(raw,barcodes)["home_results"][0]
        self.assertTrue(home["available"])
        rows=[r for r in home["receiver_half_results"] if r["receiver_culture"]=="DT2"]
        self.assertEqual([r["receiver_test_occupied_barcodes"] for r in rows],[1,1])
        raw[np.flatnonzero(parity==0),2]=0
        home=EXPLORATORY.temporal_diagnostic(raw,barcodes)["home_results"][0]
        self.assertFalse(home["available"])
        self.assertTrue(any(not r["available"] for r in home["receiver_half_results"]))

    def test_no_eligible_temporal_home_is_explicitly_unavailable(self):
        raw=invented_counts()
        result=EXPLORATORY.temporal_diagnostic(raw,np.array([f"FAKE_{i}" for i in range(len(raw))]))
        self.assertEqual(result["available_home_cultures"],0)
        self.assertIsNone(result["equal_home_mean_gain"])

    def test_highest_category_removal_renormalizes_both_predictions(self):
        obs=MODEL.prepare_observations(invented_counts())
        fit=MODEL.fit_outer_fold(obs.early,obs.late,0)
        result=ANALYSIS.top_category_removal(obs.late[0],fit)
        top=int(np.argmax(obs.late[0]));keep=[i for i in range(len(obs.late[0])) if i!=top]
        p=obs.late[0,keep];p=p/p.sum()
        q0=np.exp(fit.base_log_prediction[keep]);q0=q0/q0.sum()
        q1=np.exp(fit.full_log_prediction[keep]);q1=q1/q1.sum()
        self.assertAlmostEqual(result["gain"],math.fsum(v*math.log(f/b) for v,f,b in zip(p,q1,q0)),places=12)


if __name__ == "__main__":
    unittest.main()
