import os
from flask import Blueprint, request, redirect, url_for, render_template, current_app, session
from werkzeug.utils import secure_filename
from src.core.data_manager import get_file_columns

# Define the blueprint
data_upload_bp = Blueprint(
    'data_upload_bp',
    __name__,
    template_folder='../../templates'
)

ALLOWED_EXTENSIONS = {'csv', 'tsv', 'txt', 'xlsx'}

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@data_upload_bp.route('/upload', methods=['GET', 'POST'])
def upload_file():
    if request.method == 'POST':
        files = request.files.getlist('file')
        if not files or files[0].filename == '':
            # TODO: Add flash message for error: No files selected
            return redirect(request.url)

        upload_folder = os.path.join(current_app.root_path, '..', 'data', 'uploads')
        os.makedirs(upload_folder, exist_ok=True)

        files_data = []
        for file in files:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                file_path = os.path.join(upload_folder, filename)
                file.save(file_path)
                
                # Get columns for this specific file
                columns = get_file_columns(file_path)
                if columns:
                    files_data.append({'filename': filename, 'columns': columns})
        
        if not files_data:
            # TODO: Add flash message for error: No valid files uploaded or columns read
            return redirect(request.url)

        session['files_data'] = files_data
        
        # Redirect to the new configuration page
        return redirect(url_for('plot_api_bp.configure_plots'))
    
    return render_template('upload.html')
