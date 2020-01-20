# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin sqlcustom [app_label]'
# into your database.
from __future__ import unicode_literals

from django.db import models


class ChipData(models.Model):
    id = models.IntegerField(db_column='ID', primary_key=True)  # Field name made lowercase.
    chip_id = models.IntegerField(db_column='CHIP_ID', blank=True, null=True)  # Field name made lowercase.
    raw = models.CharField(db_column='RAW', max_length=250, blank=True, null=True)  # Field name made lowercase.
    at_100 = models.CharField(db_column='AT_100', max_length=250, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'CHIP_DATA'


class ChipSeq(models.Model):
    chip_id = models.IntegerField(db_column='CHIP_ID', primary_key=True)  # Field name made lowercase.
    tf_name = models.CharField(db_column='TF_NAME', max_length=45, blank=True, null=True)  # Field name made lowercase.
    tf = models.ForeignKey('TranscriptionFactor', db_column='TF_ID', blank=True, null=True)  # Field name made lowercase.
    database = models.CharField(db_column='DATABASE', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'CHIP_SEQ'


class Matrix(models.Model):
    id = models.IntegerField(db_column='ID', primary_key=True)  # Field name made lowercase.
    motif_id = models.CharField(db_column='MOTIF_ID', max_length=45, blank=True, null=True)  # Field name made lowercase.
    motif_name = models.CharField(db_column='MOTIF_NAME', max_length=100, blank=True, null=True)  # Field name made lowercase.
    tf = models.ForeignKey('TranscriptionFactor', db_column='TF_ID', blank=True, null=True)  # Field name made lowercase.
    collection = models.CharField(db_column='COLLECTION', max_length=45, blank=True, null=True)  # Field name made lowercase.
    link = models.ForeignKey('UrlTab', db_column='LINK_ID', blank=True, null=True)  # Field name made lowercase.
    type = models.CharField(db_column='TYPE', max_length=45, blank=True, null=True)  # Field name made lowercase.
    pub = models.ForeignKey('Publications', db_column='PUB_ID', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'MATRIX'


class MatrixData(models.Model):
    id = models.AutoField(unique=True)
    matrix = models.ForeignKey(Matrix, db_column='MATRIX_ID')  # Field name made lowercase.
    row = models.CharField(db_column='ROW', max_length=45)  # Field name made lowercase.
    col = models.IntegerField(db_column='COL')  # Field name made lowercase.
    val = models.FloatField(db_column='VAL', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'MATRIX_DATA'
        unique_together = (('id', 'MATRIX_ID', 'ROW', 'COL'),)


class Pbm(models.Model):
    pbm_id = models.IntegerField(db_column='PBM_ID', primary_key=True)  # Field name made lowercase.
    tf_name = models.CharField(db_column='TF_NAME', max_length=45, blank=True, null=True)  # Field name made lowercase.
    tf_id = models.CharField(db_column='TF_ID', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'PBM'


class PbmData(models.Model):
    id = models.IntegerField(db_column='ID', primary_key=True)  # Field name made lowercase.
    pbm = models.ForeignKey(Pbm, db_column='PBM_ID')  # Field name made lowercase.
    pbm_debru = models.CharField(db_column='PBM_DEBRU', max_length=45, blank=True, null=True)  # Field name made lowercase.
    source = models.CharField(db_column='SOURCE', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'PBM_DATA'


class Publications(models.Model):
    pub_id = models.IntegerField(db_column='PUB_ID', primary_key=True)  # Field name made lowercase.
    small_ref = models.CharField(db_column='SMALL_REF', max_length=45, blank=True, null=True)  # Field name made lowercase.
    full_ref = models.CharField(db_column='FULL_REF', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'PUBLICATIONS'


class TfClass(models.Model):
    tf_class_id = models.CharField(db_column='TF_CLASS_ID', primary_key=True, max_length=45)  # Field name made lowercase.
    tf_class = models.CharField(db_column='TF_CLASS', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'TF_CLASS'


class TranscriptionFactor(models.Model):
    tf_id = models.CharField(db_column='TF_ID', primary_key=True, max_length=45)  # Field name made lowercase.
    tf_name = models.CharField(db_column='TF_NAME', max_length=45, blank=True, null=True)  # Field name made lowercase.
    alt_tf_name = models.CharField(db_column='ALT_TF_NAME', max_length=45, blank=True, null=True)  # Field name made lowercase.
    tf_class = models.ForeignKey(TfClass, db_column='TF_CLASS_ID')  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'TRANSCRIPTION_FACTOR'


class UrlTab(models.Model):
    link_id = models.IntegerField(db_column='LINK_ID', primary_key=True)  # Field name made lowercase.
    url = models.CharField(db_column='URL', max_length=200, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'URL_TAB'


class AuthGroup(models.Model):
    name = models.CharField(unique=True, max_length=80)

    class Meta:
        managed = False
        db_table = 'auth_group'


class AuthGroupPermissions(models.Model):
    group = models.ForeignKey(AuthGroup)
    permission = models.ForeignKey('AuthPermission')

    class Meta:
        managed = False
        db_table = 'auth_group_permissions'
        unique_together = (('group_id', 'permission_id'),)


class AuthPermission(models.Model):
    name = models.CharField(max_length=255)
    content_type = models.ForeignKey('DjangoContentType')
    codename = models.CharField(max_length=100)

    class Meta:
        managed = False
        db_table = 'auth_permission'
        unique_together = (('content_type_id', 'codename'),)


class AuthUser(models.Model):
    password = models.CharField(max_length=128)
    last_login = models.DateTimeField(blank=True, null=True)
    is_superuser = models.IntegerField()
    username = models.CharField(unique=True, max_length=30)
    first_name = models.CharField(max_length=30)
    last_name = models.CharField(max_length=30)
    email = models.CharField(max_length=254)
    is_staff = models.IntegerField()
    is_active = models.IntegerField()
    date_joined = models.DateTimeField()

    class Meta:
        managed = False
        db_table = 'auth_user'


class AuthUserGroups(models.Model):
    user = models.ForeignKey(AuthUser)
    group = models.ForeignKey(AuthGroup)

    class Meta:
        managed = False
        db_table = 'auth_user_groups'
        unique_together = (('user_id', 'group_id'),)


class AuthUserUserPermissions(models.Model):
    user = models.ForeignKey(AuthUser)
    permission = models.ForeignKey(AuthPermission)

    class Meta:
        managed = False
        db_table = 'auth_user_user_permissions'
        unique_together = (('user_id', 'permission_id'),)


class DjangoAdminLog(models.Model):
    action_time = models.DateTimeField()
    object_id = models.TextField(blank=True, null=True)
    object_repr = models.CharField(max_length=200)
    action_flag = models.SmallIntegerField()
    change_message = models.TextField()
    content_type = models.ForeignKey('DjangoContentType', blank=True, null=True)
    user = models.ForeignKey(AuthUser)

    class Meta:
        managed = False
        db_table = 'django_admin_log'


class DjangoContentType(models.Model):
    app_label = models.CharField(max_length=100)
    model = models.CharField(max_length=100)

    class Meta:
        managed = False
        db_table = 'django_content_type'
        unique_together = (('app_label', 'model'),)


class DjangoMigrations(models.Model):
    app = models.CharField(max_length=255)
    name = models.CharField(max_length=255)
    applied = models.DateTimeField()

    class Meta:
        managed = False
        db_table = 'django_migrations'


class DjangoSession(models.Model):
    session_key = models.CharField(primary_key=True, max_length=40)
    session_data = models.TextField()
    expire_date = models.DateTimeField()

    class Meta:
        managed = False
        db_table = 'django_session'
