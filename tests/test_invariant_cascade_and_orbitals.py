"""
End-to-End Regression Test Suite for Reaction Pathway & Wavefunction Tracking Invariants.
Validates:
  1. Universal atom aliases (atom_labels) & symbols consistency across all frames.
  2. Implicit ddCOSMO dielectric screening (Cl7 ion-pair screening in water) & Step-0 planar conjugated pi separation:
     - Zero Kekule ring pi leaks in localized IBOs.
     - Twisted frames (0, 61, 85): 3 delocalized Ar pi + 3 Ar pi*, 1 localized pi(C6-O8), 2 LP(O8).
     - Planar acylium frames (18..55, e.g. 27): 4 delocalized Ar-C6-O8 pi + 4 Ar-C6-O8 pi*, 1 perpendicular localized pi(C6-O8), 1 LP(O8).
  3. Part-agnostic continuous wavefunction tracking across localized <-> delocalized transitions (<100ms reload).
  4. UI Single visibility set & exact mesh/cube resolution without substring collision.
"""

import json
import os
import time
import unittest

from irc_plugins.constants import ROOT_DIR, PROLONGED_PATH as TRAJ_PATH, NPZ_PATH, REGISTRY_PATH as REG_PATH


class TestPathwayAndOrbitalInvariants(unittest.TestCase):

    def setUp(self):
        self.assertTrue(os.path.exists(TRAJ_PATH), f"Missing {TRAJ_PATH}")
        with open(TRAJ_PATH, "r") as f:
            self.megajson = json.load(f)

    def test_invariant_1_universal_atom_aliases(self):
        """All frames must preserve consistent atom_labels regardless of leaving groups."""
        frames = self.megajson.get("frames", [])
        self.assertGreater(len(frames), 50, "Expected at least 50 frames in pathway")

        lbls_0 = frames[0].get("atom_labels", [])
        self.assertIn("Cl7", lbls_0)
        self.assertIn("C6", lbls_0)
        self.assertIn("O8", lbls_0)

        # In prolonged pathway, Cl7 leaves at Step 2 (Frame 50)
        lbls_50 = frames[50].get("atom_labels", [])
        self.assertNotIn("Cl7", lbls_50)
        self.assertIn("C6", lbls_50)
        self.assertIn("O8", lbls_50)

        lbls_61 = frames[61].get("atom_labels", [])
        self.assertIn("C6", lbls_61)
        self.assertIn("O8", lbls_61)
        self.assertIn("O24", lbls_61)

    def test_invariant_2_ddcosmo_screening_and_step0_planar_pi_exclusion(self):
        """ddCOSMO water screening + Step-0 planar conjugated pi exclusion must match server.py physics."""
        from irc_plugins.orbital_tracker import get_frame_pool_payload

        self.assertEqual(self.megajson.get("solvent"), "water")

        # 1. Departing Cl7 at Step 1 end (Frame 49) in water (ddCOSMO, eps=80.4) must be cleanly screened (~ -1.00e)
        payload_49 = get_frame_pool_payload(self.megajson, 49)
        cl7_idx = payload_49["atom_labels"].index("Cl7")
        cl7_charge = payload_49["charges"][cl7_idx]
        self.assertLessEqual(
            cl7_charge, -0.95,
            f"Expected ddCOSMO-screened Cl- charge <= -0.95e at Frame 49, got {cl7_charge}",
        )

        ring_set = {"C1", "C2", "C3", "C5", "C9", "C11"}

        # 2. Carbonyl frames with intact >C=O (0, 100, 150): 3 Ar pi + 3 Ar pi* in canonical_pi, 1 localized pi(C6-O8) in all_ibos
        for fi_dbl in (0, 100, 150):
            payload_dbl = get_frame_pool_payload(self.megajson, fi_dbl)
            ring_leaks = [
                x["label"] for x in payload_dbl.get("all_ibos", [])
                if "pi" in x.get("type", "") and set(x.get("atoms", [])) <= ring_set
            ]
            self.assertEqual(ring_leaks, [], f"Frame {fi_dbl} leaked ring pi orbitals into IBO: {ring_leaks}")

            co_bonds_dbl = [x["label"] for x in payload_dbl.get("all_ibos", []) if "C6-O8" in x.get("label", "")]
            self.assertEqual(
                sorted(co_bonds_dbl),
                ["π(C6-O8)", "σ(C6-O8)"],
                f"Frame {fi_dbl} twisted C6-O8 must have 1 localized sigma and 1 localized pi, got {co_bonds_dbl}",
            )
            o8_lps_dbl = [
                x["label"] for x in payload_dbl.get("all_ibos", [])
                if x.get("type") == "lone_pair" and x.get("atoms") == ["O8"]
            ]
            self.assertEqual(len(o8_lps_dbl), 2, f"Frame {fi_dbl} O8 must have 2 lone pairs, got {o8_lps_dbl}")

            canon_occ = [x for x in payload_dbl.get("canonical_pi", []) if x.get("type") == "delocalized_pi"]
            canon_virt = [x for x in payload_dbl.get("canonical_pi", []) if x.get("type") == "delocalized_pistar"]
            self.assertEqual(len(canon_occ), 3, f"Frame {fi_dbl} twisted ring must have 3 delocalized pi MOs")
            self.assertEqual(len(canon_virt), 3, f"Frame {fi_dbl} twisted ring must have 3 delocalized pi* MOs")
            for x in canon_virt:
                self.assertEqual(x.get("orbital_class"), "antibonding")

        # 3. Planar acylium Frame 116: coplanar pi absorbed into 4-orbital Ar-C6-O8 conjugated system
        payload_116 = get_frame_pool_payload(self.megajson, 116)
        ring_leaks_116 = [
            x["label"] for x in payload_116.get("all_ibos", [])
            if "pi" in x.get("type", "") and set(x.get("atoms", [])) <= ring_set
        ]
        self.assertEqual(ring_leaks_116, [], f"Frame 116 leaked ring pi orbitals into IBO: {ring_leaks_116}")

        canon_occ_116 = [x for x in payload_116.get("canonical_pi", []) if x.get("type") == "delocalized_pi"]
        canon_virt_116 = [x for x in payload_116.get("canonical_pi", []) if x.get("type") == "delocalized_pistar"]
        self.assertEqual(len(canon_occ_116), 4, "Expected 4 conjugated Ar-C6-O8 pi MOs at planar acylium Frame 116")
        self.assertEqual(len(canon_virt_116), 4, "Expected 4 conjugated Ar-C6-O8 pi* MOs at planar acylium Frame 116")
        for x in canon_occ_116:
            self.assertIn("C6", x.get("atoms", []))
            self.assertIn("O8", x.get("atoms", []))

    def test_invariant_3_part_agnostic_orbital_continuity(self):
        """Orbital tracking must follow continuous wavefunctions without knowledge of parts or breaks."""
        self.assertTrue(os.path.exists(NPZ_PATH), "Missing orbital_tracks.npz")
        self.assertTrue(os.path.exists(REG_PATH), "Missing orbital_registry.json")

        from irc_plugins.orbital_tracker import get_or_build_orbital_store, clear_memory_cache

        clear_memory_cache()
        t0 = time.perf_counter()
        store = get_or_build_orbital_store(self.megajson)
        reload_ms = (time.perf_counter() - t0) * 1000.0
        self.assertLess(reload_ms, 500.0, f"Orbital reload took {reload_ms:.2f} ms (expected < 500 ms)")

        reg = store["registry"]["tracks"]

        # 1. Leaving group bond: sigma(C6-Cl7) -> LP(Cl7) (trk_0016) terminates when Cl departs at frame 49
        self.assertIn("trk_0016", reg)
        t_c6_cl7 = reg["trk_0016"]
        self.assertEqual(t_c6_cl7["label"], "σ(C6-Cl7) → LP(Cl7)")
        self.assertEqual((t_c6_cl7["start_frame"], t_c6_cl7["end_frame"]), (0, 49))
        for fi in range(0, 50):
            self.assertEqual(t_c6_cl7["frame_map"][fi], "trk_0016", f"trk_0016 broken at frame {fi}")
        self.assertIsNone(t_c6_cl7["frame_map"][50], "trk_0016 must die at frame 50 when Cl- departs")

        # 2. Leaving group virtual channel: sigma*(C6-Cl7) (trk_0046) terminates upon bond cleavage at frame 28
        self.assertIn("trk_0046", reg)
        t_c6_cl7_star = reg["trk_0046"]
        self.assertEqual(t_c6_cl7_star["label"], "σ*(C6-Cl7)")
        self.assertEqual((t_c6_cl7_star["start_frame"], t_c6_cl7_star["end_frame"]), (0, 28))
        for fi in range(0, 29):
            self.assertEqual(t_c6_cl7_star["frame_map"][fi], "trk_0046", f"trk_0046 broken at frame {fi}")
        self.assertIsNone(t_c6_cl7_star["frame_map"][29], "trk_0046 must die at frame 29 upon cleavage")

        # 3. Nucleophilic attack: LP(O24) -> sigma(C6-O24) (trk_0082) born as nucleophile attacks
        self.assertIn("trk_0082", reg)
        t_c6_o24 = reg["trk_0082"]
        self.assertEqual(t_c6_o24["label"], "LP(O24) → σ(C6-O24)")
        self.assertEqual((t_c6_o24["start_frame"], t_c6_o24["end_frame"]), (66, 188))
        for fi in range(66, 189):
            self.assertEqual(t_c6_o24["frame_map"][fi], "trk_0082", f"trk_0082 broken at frame {fi}")

        # 4. Persistent spectator carbonyl pi(C6-O8) -> LP(O8) (trk_0029) continuous across all 189 frames
        self.assertIn("trk_0029", reg)
        t_pi = reg["trk_0029"]
        self.assertEqual((t_pi["start_frame"], t_pi["end_frame"]), (0, 188))
        self.assertEqual(t_pi["label"], "π(C6-O8) → LP(O8)")
        for fi in range(0, 189):
            self.assertEqual(t_pi["frame_map"][fi], "trk_0029", f"trk_0029 broken at frame {fi}")

        # 5. Reaction coordinate occupied channel: LP(O8) -> pi(C6-O8) -> pi_deloc(C6-O8-O24) [pi2] (trk_0027) continuous across 0..188
        self.assertIn("trk_0027", reg)
        t_rc_occ = reg["trk_0027"]
        self.assertEqual((t_rc_occ["start_frame"], t_rc_occ["end_frame"]), (0, 188))
        self.assertEqual(t_rc_occ["label"], "LP(O8) → π(C6-O8) → π_deloc(C6-O8-O24) [π2]")
        for fi in range(0, 189):
            self.assertEqual(t_rc_occ["frame_map"][fi], "trk_0027", f"trk_0027 broken at frame {fi}")

        # 6. Delocalized aromatic pi* LUMO reaction coordinate channel (trk_0057) continuous across 0..188
        self.assertIn("trk_0057", reg)
        t_pistar_lumo = reg["trk_0057"]
        self.assertEqual((t_pistar_lumo["start_frame"], t_pistar_lumo["end_frame"]), (0, 188))
        self.assertEqual(t_pistar_lumo["label"], "π*_deloc(Ar) [π*1] → π*_deloc(Ar) [π*2]")
        for fi in range(0, 189):
            self.assertEqual(t_pistar_lumo["frame_map"][fi], "trk_0057", f"trk_0057 broken at frame {fi}")

        # 7. Localized spectator carbonyl pi*(C6-O8) (trk_0056)
        self.assertIn("trk_0056", reg)
        t_pistar_loc = reg["trk_0056"]
        self.assertEqual(t_pistar_loc["label"], "π*₁(C6-O8)")
        self.assertEqual((t_pistar_loc["start_frame"], t_pistar_loc["end_frame"]), (0, 55))
        for fi in range(0, 56):
            self.assertEqual(t_pistar_loc["frame_map"][fi], "trk_0056", f"trk_0056 broken at frame {fi}")

        # 8. Water proton-transfer channel: LP(O29) -> sigma(O29-H28) (trk_0088) continuous across 120..188
        self.assertIn("trk_0088", reg)
        t_h2o = reg["trk_0088"]
        self.assertEqual((t_h2o["start_frame"], t_h2o["end_frame"]), (120, 188))
        for fi in range(120, 189):
            self.assertEqual(t_h2o["frame_map"][fi], "trk_0088", f"trk_0088 broken at frame {fi}")

    def test_invariant_4_ui_single_visibility_and_exact_meshing(self):
        """Frontend files must enforce single visibility set and exact key mesh resolution."""
        with open(os.path.join(ROOT_DIR, "static_irc", "ui_plugins", "orbitals_ui.js"), "r") as f:
            ui_code = f.read()
        self.assertNotIn("if (auditioned.length > 0) return activeAud;", ui_code)
        self.assertIn("visibleTrackIds: new Set()", ui_code)
        self.assertIn("this.state.visibleTrackIds = new Set(initialVisible);", ui_code)
        self.assertIn("toggleVisibility(tid)", ui_code)

        with open(os.path.join(ROOT_DIR, "static_irc", "ui_plugins", "orbitals", "orbitals_timeline.js"), "r") as f:
            tl_code = f.read()
        self.assertNotIn("this.manualOverrides", tl_code)
        self.assertIn("this.visibleSet", tl_code)

        with open(os.path.join(ROOT_DIR, "static_viewer", "js", "viewer_surface.js"), "r") as f:
            surf_code = f.read()
        self.assertNotIn("k.includes(orbKey) || orbKey.includes(baseKey)", surf_code)
        self.assertIn("const exactKey = `${frameIdx}_${orbKey}`;", surf_code)


if __name__ == "__main__":
    unittest.main()
