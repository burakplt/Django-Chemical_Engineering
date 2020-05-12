from django.urls import path, re_path
from . import views


urlpatterns = [
    path('', views.vle_home, name="vle_home"),
    #path('bubble-dew', views.bub_dew, name="bub_dew"),
    path('fugacity', views.fugacity, name="fugacity"),
    path('fugacity/activity', views.activity, name="activity"),
    path('pt-flash', views.PT_flash, name="pt_flash"),
    #re_path(r'^ajax/bubble-dew-result/$', views.bub_dew_result, name="bubble_dew_result"),
    re_path(r'^ajax/fugacity-result/$', views.fugacity_vapor_result, name="fugacity_vapor_result"),
    re_path(r'^ajax/activity-result/$', views.activity_result, name="activity_result"),
    
] 