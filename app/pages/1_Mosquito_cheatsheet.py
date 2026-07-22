from __future__ import annotations

from pathlib import Path

import streamlit as st


ASSET_DIR = Path(__file__).resolve().parents[1] / "assets"


st.set_page_config(
    page_title="Mosquito basics · Aedes RNA Atlas",
    page_icon="🦟",
    layout="wide",
    initial_sidebar_state="collapsed",
)

st.markdown(
    """
    <style>
    .block-container { padding-top: 2.4rem; padding-bottom: 3rem; max-width: 1180px; }
    .cheat-lead { color: #9aa8a5; max-width: 850px; margin-bottom: 1rem; }
    </style>
    """,
    unsafe_allow_html=True,
)

st.page_link("app.py", label="Back to expression explorer", icon="⬅️")
st.title("Mosquito basics")
st.markdown(
    '<div class="cheat-lead">Plain-language definitions for the tissues and adult reproductive conditions used in the RNA Atlas.</div>',
    unsafe_allow_html=True,
)

st.info(
    "Both RNA-seq studies in this atlas sampled adult mosquitoes. Egg, larva, and pupa are included below for orientation, but they are not expression conditions in these matrices."
)

st.image(
    ASSET_DIR / "mosquito_body_parts_reference.png",
    caption=(
        "General mosquito anatomy reference. The source image is labeled Culex pipiens, "
        "so use it for broad body-part orientation rather than Aedes-specific morphology."
    ),
    width=600,
)

st.markdown("## How to read a plot row")
st.markdown(
    "**`Antenna · Blood-fed`** means RNA was measured from antenna tissue dissected from blood-fed adult females. "
    "The word before `·` is the tissue; the phrase after it is the sex, feeding state, or reproductive time point."
)

st.markdown("## Body parts used in the atlas")
st.caption("🧪 marks a tissue represented directly in the RNA-seq matrices. Other parts are included so the anatomy makes sense.")
head, thorax, abdomen = st.columns(3)

with head:
    with st.container(border=True):
        st.markdown("### Head · sensing and feeding")
        st.markdown(
            "**🧪 Antenna / antennae**  \n"
            "Paired sensory appendages that detect odors and air cues."
        )
        st.markdown(
            "**Compound eyes**  \n"
            "The two large eyes on the head; built from many visual units and specialized for detecting movement."
        )
        st.markdown(
            "**🧪 Brain**  \n"
            "Central nervous tissue that receives and processes sensory information."
        )
        st.markdown(
            "**🧪 Palps / maxillary palps**  \n"
            "Short paired sensory appendages beside the proboscis; important for odor sensing."
        )
        st.markdown(
            "**🧪 Proboscis**  \n"
            "The elongated feeding mouthpart. Females use it to take blood meals; both sexes feed on sugars."
        )
        st.markdown(
            "**🧪 Rostrum**  \n"
            "In the 2016 study, this was a combined dissection containing the maxillary palps and proboscis."
        )

with thorax:
    with st.container(border=True):
        st.markdown("### Thorax · movement")
        st.markdown(
            "**Wings**  \n"
            "One pair of wings used for flight."
        )
        st.markdown(
            "**Halteres**  \n"
            "Small balancing organs behind the wings that help stabilize flight."
        )
        st.markdown(
            "**🧪 Forelegs**  \n"
            "The front pair of legs."
        )
        st.markdown(
            "**🧪 Midlegs**  \n"
            "The middle pair of legs."
        )
        st.markdown(
            "**🧪 Hindlegs**  \n"
            "The rear pair of legs."
        )
        st.markdown(
            "**Leg segments**  \n"
            "From body outward, the prominent sections include the **femur**, **tibia**, and multi-part **tarsus** (the foot-like end). The atlas samples whole leg pairs, not individual segments."
        )
        st.caption(
            "All six legs and the wings attach to the thorax. Wings and whole thorax were not separate samples in these matrices."
        )

with abdomen:
    with st.container(border=True):
        st.markdown("### Abdomen · digestion and reproduction")
        st.markdown(
            "**Abdominal segments**  \n"
            "The abdomen is the long, visibly segmented rear body region; it contains much of the digestive and reproductive system."
        )
        st.markdown(
            "**🧪 Ovaries**  \n"
            "Female reproductive organs in which eggs develop and can be retained. The drought study measured ovaries only."
        )
        st.markdown(
            "**🧪 Abdominal tip**  \n"
            "The final three abdominal segments. The study dissection included genitalia and the female ovipositor."
        )
        st.markdown(
            "**Ovipositor**  \n"
            "The egg-laying structure at the end of the female abdomen. It is part of the abdominal-tip sample, not a separate tissue row."
        )
        st.markdown(
            "**Genitalia**  \n"
            "Reproductive structures at the tip of the abdomen; included inside the abdominal-tip dissection."
        )

st.markdown("## Mosquito life stages")
egg, larva, pupa, adult = st.columns(4)
with egg:
    with st.container(border=True):
        st.markdown("### 1 · Egg")
        st.write("Laid near water. Aedes eggs can tolerate drying and hatch after being covered by water.")
with larva:
    with st.container(border=True):
        st.markdown("### 2 · Larva")
        st.write("Aquatic, feeding, growing stage—often called a wiggler.")
with pupa:
    with st.container(border=True):
        st.markdown("### 3 · Pupa")
        st.write("Aquatic, non-feeding transition stage from larva to adult.")
with adult:
    with st.container(border=True):
        st.markdown("### 4 · Adult")
        st.write("Flying stage. This is the only life stage represented in the current RNA Atlas datasets.")

st.markdown("## Adult feeding and reproductive states")
st.markdown(
    """
| Atlas label | Plain-language meaning in the 2016 neurotranscriptome study |
| --- | --- |
| **Non-blood-fed / sugar-fed** | Baseline adult females that had not taken a blood meal and could feed on sugar. Sample code: `SF`. |
| **Blood-fed** | Adult females sampled 48 hours after a blood meal. Sample code: `BF`. |
| **Gravid / oviposition-stage** | Adult females sampled 96 hours after a blood meal, carrying developed eggs and near the egg-laying stage. Sample code: `O`. |
| **Male** | Adult male tissue; this is a sex comparison, not a female reproductive state. |
"""
)

st.markdown("## Drought-study ovary timeline")
st.markdown(
    "The drought-resilience study follows adult female ovaries before and after a blood meal. "
    "Every label is elapsed time after that meal; **PBM** means **post-blood-meal**."
)
st.code(
    "Non-blood-fed → 3 h → 6 h → 12 h → 24 h → 48 h → 72 h → 96 h\n"
    "                                      ↘ day 6: eggs retained\n"
    "                                      ↘ day 6: eggs laid <5 h prior\n"
    "                                      ↘ day 13: eggs laid >1 week prior",
    language=None,
)
st.markdown(
    """
| Special label | What it means |
| --- | --- |
| **6 days PBM · eggs retained** | Freshwater for egg laying was unavailable, so mature eggs remained in the ovaries. This models a drought-like condition. |
| **6 days PBM · eggs laid <5 hours prior** | Females had very recently laid their eggs; ovaries were sampled within five hours afterward. |
| **13 days PBM · eggs laid >1 week prior** | A later post-egg-laying ovary state, more than a week after laying. |
"""
)

st.markdown("## Quick vocabulary")
st.markdown(
    """
- **Gravid:** carrying developed eggs.
- **Oviposition:** laying eggs.
- **Post-blood-meal / PBM:** elapsed time since the female took a blood meal.
- **Replicate:** an independently prepared biological sample; the small plot points are replicates.
- **Group median:** the middle expression value among replicates; shown as an orange diamond.
- **TPM:** transcripts per million, a within-sample normalized measure of transcript abundance.
"""
)

st.markdown("## Sources")
st.markdown(
    "- [Matthews et al. 2016 · neurotranscriptome study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4704297/)  \n"
    "- [Venkataraman et al. 2023 · drought-resilience study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10076016/)  \n"
    "- [CDC · parts of an adult mosquito](https://www.cdc.gov/mosquitoes/about/index.html)  \n"
    "- [CDC · life cycle of Aedes mosquitoes](https://www.cdc.gov/mosquitoes/about/life-cycle-of-aedes-mosquitoes.html)"
)
