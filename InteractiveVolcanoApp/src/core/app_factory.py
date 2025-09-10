from flask import Flask, render_template, session, current_app
import os
import json

def create_app():
    """Create and configure an instance of the Flask application."""
    app = Flask(__name__, template_folder='../../templates', static_folder='../../static')

    # Set a secret key for session management.
    app.config['SECRET_KEY'] = os.urandom(24)

    # Register blueprints
    from src.api.data_upload import data_upload_bp
    from src.api.plot_api import plot_api_bp
    from src.api.subtiwiki_api_bp import subtiwiki_api_bp
    from src.api.selection_api_bp import selection_api_bp
    app.register_blueprint(data_upload_bp)
    app.register_blueprint(plot_api_bp)
    app.register_blueprint(subtiwiki_api_bp, url_prefix='/api')
    app.register_blueprint(selection_api_bp, url_prefix='/api')

    def get_plot_data():
        """Helper function to read plot data from a temporary file.

        Robust version: keeps session key and does not delete the file so
        refreshes can still load the same data during the session.
        """
        plot_data = []
        # Use get() to persist the key across page reloads
        processed_file_id = session.get('processed_data_file', None)
        if processed_file_id:
            processed_folder = os.path.join(current_app.root_path, '..', 'data', 'processed')
            file_path = os.path.join(processed_folder, processed_file_id)
            try:
                with open(file_path, 'r') as f:
                    plot_data = json.load(f)
                # Do NOT delete here; let a periodic cleanup handle stale files
                # os.remove(file_path)
            except Exception as e:
                print(f"Error reading processed data file: {e}")
        return plot_data

    @app.route('/')
    def index():
        """Render the main analysis page, displaying plots if data is available."""
        plot_data = get_plot_data()
        return render_template('analysis.html', plot_data=plot_data)

    @app.route('/analysis')
    def analysis():
        """Render the main analysis page, displaying plots if data is available."""
        plot_data = get_plot_data()
        return render_template('analysis.html', plot_data=plot_data)

    return app
