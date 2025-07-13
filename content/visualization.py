"""
Visualization Page

This page will contain various data visualizations for peptide analysis
"""

import streamlit as st
from src.common.common import page_setup, v_space


def main():
    """Main function for the visualization page."""
    page_setup()

    st.title("📊 Visualization")

    # Under development notice
    st.info("🚧 This page is currently under development")


# Execute the main function
main()
