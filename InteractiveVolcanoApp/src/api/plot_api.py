import os
import uuid
import json
from collections import defaultdict
from flask import Blueprint, render_template, session, redirect, url_for, request, current_app
from src.core.data_manager import create_plot_data

# Define the blueprint
plot_api_bp = Blueprint(
    'plot_api_bp',
    __name__,
    template_folder='../../templates'
)

@plot_api_bp.route('/configure', methods=['GET'])
def configure_plots():
    files_data = session.get('files_data', [])
    if not files_data:
        return redirect(url_for('data_upload_bp.upload_file'))
    return render_template('configure.html', files_data=files_data)

@plot_api_bp.route('/process', methods=['POST'])
def process_plots():
    print("=== PROCESS PLOTS ENDPOINT CALLED ===")
    form_data = request.form
    print(f"Form data received: {dict(form_data)}")
    files_data = session.get('files_data', [])
    print(f"Files data from session: {files_data}")

    # Safety check: if session is empty, redirect to upload
    if not files_data:
        return redirect(url_for('data_upload_bp.upload_file'))

    parsed_configs = defaultdict(lambda: defaultdict(dict))
    for key, value in form_data.items():
        parts = key.split('_')
        file_index = int(parts[1]) - 1
        plot_index = int(parts[3]) - 1
        col_type_key = '_'.join(parts[4:])
        parsed_configs[file_index][plot_index][col_type_key] = value
        print(f"🔧 Parsed config: file_{file_index}, plot_{plot_index}, {col_type_key} = '{value}'")

    plot_data_list = []
    upload_folder = os.path.join(current_app.root_path, '..', 'data', 'uploads')

    for file_index in sorted(parsed_configs.keys()):
        filename = files_data[file_index]['filename']
        for plot_index in sorted(parsed_configs[file_index].keys()):
            config = parsed_configs[file_index][plot_index]
            plot_data = create_plot_data(filename, config, upload_folder)
            if plot_data:
                plot_data_list.append(plot_data)

    print(f"Created {len(plot_data_list)} plot data entries")
    
    if plot_data_list:
        # Save the processed data to a temporary file instead of the session
        processed_data_id = f"plot_data_{uuid.uuid4()}.json"
        processed_folder = os.path.join(current_app.root_path, '..', 'data', 'processed')
        os.makedirs(processed_folder, exist_ok=True)
        processed_file_path = os.path.join(processed_folder, processed_data_id)

        with open(processed_file_path, 'w') as f:
            json.dump(plot_data_list, f)

        # Store only the ID in the session
        session['processed_data_file'] = processed_data_id
        print(f"Saved plot data to: {processed_data_id}")
    else:
        print("WARNING: No plot data was created!")

    print("Redirecting to analysis page...")
    return redirect(url_for('analysis'))
