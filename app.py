import streamlit as st
import matplotlib.pyplot as plt
from pagina.page_primers import page_primer

# --------- Configure Page --------
st.set_page_config(page_title='SmartHTPCloning', layout='wide')

with open('style.css') as f:
    st.markdown(f'<style>{f.read()}<\style>', unsafe_allow_html=True)

# ------------- Menu ------------
st.markdown("<h1 style = 'text-align:center'> SmartHTPCloning </h3>", unsafe_allow_html=True)

# --------- Configure DNA plot -------------

plt.rcParams.update({
    "figure.facecolor": (1.0, 1.0, 1.0, 0.0),  # red   with alpha = 30%
    "axes.facecolor": (1.0, 1.0, 1.0, 0.0),  # green with alpha = 50%
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),  # blue  with alpha = 20%
})

# --- Attributes ---
page_primer()
